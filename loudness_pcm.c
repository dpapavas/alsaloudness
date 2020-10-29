/* alsaloundess, a loudness compensation plugin for ALSA.
 * Copyright (C) 2018 Papavasileiou Dimitris
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#define _POSIX_C_SOURCE 200809L

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <complex.h>
#include <pthread.h>
#include <sys/types.h>

#ifndef NDEBUG
#include <time.h>
#endif

#include <fftw3.h>

#include <alsa/asoundlib.h>
#include <alsa/pcm.h>
#include <alsa/pcm_external.h>
#include <alsa/error.h>
#include <alsa/control.h>
#include <alsa/sound/tlv.h>
#include <linux/soundcard.h>

#include "contours.h"

#ifndef NDEBUG
#define TRACE(format, ...) fprintf(stderr, "%s:%d:(%s) " format, \
                                   __FILE__, __LINE__, __func__, ##__VA_ARGS__)
#else
#define TRACE(format, ...)
#endif

#define ATTENUATION 0
#define REFERENCE 1
#define COMPENSATE 2

static const struct{
    const char *name;
    int range[2];
} CONTROLS[] = {
    [ATTENUATION] = {"Attenuation Playback Volume", {-50, 0}},
    [REFERENCE] = {"Reference Playback Volume", {50, 90}},
    [COMPENSATE] = {"Attenuation Playback Switch", {0, 1}},
};

#define CONTROLS_COUNT (sizeof(CONTROLS) / sizeof(CONTROLS[0]))

struct context {
    snd_pcm_extplug_t ext;
    snd_ctl_t *ctl;

    char *prefix;

    unsigned int impulse_length;
    unsigned int fft_length;
    unsigned int interpolation_length;

    /* Filtering-related stuff. */

    float *input, *output;
    fftwf_complex *bins;

    fftwf_plan input_to_bins, bins_to_output;

    /* Design-related stuff. */

    float *kernel;
    fftwf_complex *spline, *weights;

    fftwf_plan spline_to_kernel, kernel_to_weights;

    int values[CONTROLS_COUNT];

    const char* wisdom_path;

#ifdef WITH_THREADS
    int threads;
#endif

#ifndef NDEBUG
    unsigned long benchmarks[2];
#endif
};

static float I_0(float x)
{
    const float x_2 = x / 2;
    float t_i, s_i, s_i_minus_1, y;
    int i;

    /* Calculate the modified Bessel function I_0, by summing terms of
     * the form [(x / 2)^2]^i / (i^2)!. */

    for (i = 1, t_i = 1, s_i_minus_1 = 0, s_i = 1;
         s_i != s_i_minus_1;
         s_i_minus_1 = s_i, y = x_2 / i, i += 1, t_i *= y * y, s_i += t_i);

    return s_i;
}

static float kaiser(float i, float l, float beta)
{
    const float x = 2 * i / (l - 1) - 1;

    return I_0(beta * sqrt(1 - x * x)) / I_0(beta);
}

static void update_filter_weights(struct context *context)
{
    const unsigned int M = context->impulse_length;
    const unsigned int L = context->interpolation_length;
    const float delta_f = (float)context->ext.rate / L;

    unsigned int i, j;
    float f;

#ifndef NDEBUG
    struct timespec t_0;

    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &t_0);
#endif

    /* Sample the difference of the reference and target curve splines
     * uniformly. */

    for (i = 0, j = 0, f = 0 ; i < L / 2 + 1 ; i += 1, f += delta_f) {
        const int l_0 = context->values[REFERENCE];
        const int delta_l = context->values[ATTENUATION];

        float x, x_2, x_3, *C_ref, *C, log_f;

        /* Splines are parameterized on log(frequency), in order to
         * achieve uniform sampling across octaves and avoid vanishing
         * coefficients. */

        log_f = logf(f == 0 ? 0.1 : f);

        /* Figure out which segment to sample. */

        while (log_f > contours.knots[j + 1]) {
            j += 1;
        }

        /* Evaluate the cubic polynomial. */

        x = log_f - contours.knots[j];
        x_2 = x * x;
        x_3 = x_2 * x;

        C_ref = contours.coefficients[l_0][j];
        C = contours.coefficients[l_0 + delta_l][j];

        context->spline[i] = expf(0.05 * logf(10) *
                                  ((C[3] - C_ref[3]) +
                                   (C[2] - C_ref[2]) * x +
                                   (C[1] - C_ref[1]) * x_2 +
                                   (C[0] - C_ref[0]) * x_3));
    }

#ifdef PLOT_FILTER
    {
        FILE *fp;

        fp = fopen("sampled", "w");

        for (i = 0, j = 0, f = 0 ; i < L / 2 + 1 ; i += 1, f += delta_f) {
            fprintf (fp, "%f %f\n",
                     f == 0 ? 0.1 : f, creal(context->spline[i]));
        }

        fclose(fp);
}
#endif

    /* Transform the sampled response to the time domain and shift,
     * truncate and window it, to yield a centered impulse response of
     * the desired length. */

    fftwf_execute(context->spline_to_kernel);

    {
        float l = log2(M);
        const float beta = ((7.5 * l * l - 190.5 * l + 1296) *
                            (context->values[ATTENUATION] * -0.05 + 2.5) / 112);

        for (i = 0 ; i < M / 2; i += 1) {
            context->kernel[M / 2 + i] =
                context->kernel[i] / L * kaiser(i + M / 2, M, beta);
        }

        for (i = 0 ; i < M / 2; i += 1) {
            context->kernel[i] =
                context->kernel[L - M / 2 + i] / L * kaiser(i, M, beta);
        }

#ifdef PLOT_FILTER
        {
            FILE *fp;

            fp = fopen("parameters", "w");

            fprintf(fp, "%u %u %u %u %d %d %f\n",
                    context->impulse_length, context->fft_length,
                    context->interpolation_length, context->ext.rate,
                    context->values[REFERENCE], context->values[ATTENUATION],
                    beta);

            fclose(fp);
        }
#endif
    }

    memset(context->kernel + M, 0, (L - M) * sizeof(context->kernel[0]));

    /* Transform the truncated impulse response, zero-padded to the FFT
     * length, into the frequency domain, to yield the frequency
     * response. */

    fftwf_execute(context->kernel_to_weights);

#ifdef PLOT_FILTER
    {
        FILE *fp;
        const unsigned int N = context->fft_length;
        const float delta_f = (float)context->ext.rate / N;

        fp = fopen("realized", "w");

        for (i = 0, j = 0, f = 0 ; i < N / 2 + 1 ; i += 1, f += delta_f) {
            fprintf (fp, "%f %f %f\n",
                     f,
                     crealf(context->weights[i]),
                     cimagf(context->weights[i]));
        }

        fclose(fp);

        fp = fopen("filter.gnu", "w");

        fprintf(fp,
                "set grid mxtics xtics ytics\n"
                "set logscale x\n"
                "set xrange [10:%d]\n"
                "set key bottom left\n"
                "plot \"realized\" using 1:(20*log10(abs($2+$3*{0,1})))\\\n"
                "                  with lines\\\n"
                "                  linecolor \"#E64501\"\\\n"
                "                  title \"Realized response\",\\\n"
                "     \"sampled\" using 1:(20*log10($2))\\\n"
                "                      with lines\\\n"
                "                      dashtype \".\"\\\n"
                "                      linecolor \"#88AD41\"\\\n"
                "                      title \"Desired response\"\n",
                context->ext.rate / 2);

        fclose(fp);
    }
#endif

#ifndef NDEBUG
    {
        struct timespec t;

        clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &t);
        TRACE("(Re)calculated filter response (impulse length %d),"
              " in %.2f ms.\n",
              M,
              (t.tv_sec - t_0.tv_sec) * 1e3 +
              (t.tv_nsec - t_0.tv_nsec) * 1e-6);
    }
#endif
}

static snd_pcm_sframes_t transfer_callback(
    snd_pcm_extplug_t *ext,
    const snd_pcm_channel_area_t *dst_areas,
    snd_pcm_uframes_t dst_offset,
    const snd_pcm_channel_area_t *src_areas,
    snd_pcm_uframes_t src_offset,
    snd_pcm_uframes_t size)
{
    struct context *context = (struct context *)ext->private_data;
    float c;
    unsigned int i;
    int j;

    const unsigned int M = context->impulse_length;
    const unsigned int N = context->fft_length;
    const int C = ext->channels;

    /* Handle any pending ctl events and update the filter, if
     * necessary. */

    if (context->ctl) {
        snd_ctl_event_t *event;
        int update = 0;

        snd_ctl_event_alloca(&event);

        while (snd_ctl_read(context->ctl, event) > 0) {
            unsigned int i;

            snd_ctl_read(context->ctl, event);

            /* There's only this one type of event, at the time of
             * writing, but one never knows what the future holds. */

            if (snd_ctl_event_get_type(event) != SND_CTL_EVENT_ELEM) {
                continue;
            }

            /* Loop through our controls and see if the event concerns us.
             * Update the filter parameters, if necessary. */

            for (i = 0, update = 0 ; i < CONTROLS_COUNT ; i += 1) {
                int p;

                if (!context->prefix) {
                    p = !strcmp(CONTROLS[i].name,
                                snd_ctl_event_elem_get_name(event));
                } else {
                    char s[strlen(context->prefix) +
                           strlen(CONTROLS[i].name) + 1];

                    strcpy(s, context->prefix);
                    strcat(s, CONTROLS[i].name);

                    p = !strcmp(s, snd_ctl_event_elem_get_name(event));
                }

                if (p) {
                    snd_ctl_elem_id_t *id;
                    snd_ctl_elem_value_t *v;
                    long l;

                    snd_ctl_elem_value_alloca(&v);
                    snd_ctl_elem_id_alloca(&id);

                    snd_ctl_event_elem_get_id(event, id);
                    snd_ctl_elem_value_set_id(v, id);

                    if (snd_ctl_elem_read(context->ctl, v) == 0 &&
                        ((l = snd_ctl_elem_value_get_integer(v, 0)) !=
                         context->values[i])) {
                        context->values[i] = l;
                        update = 1;
                    }
                }
            }
        }

        /* If the filter's configuration changed, update its frequency
         * response. */

        if (update &&
            context->values[COMPENSATE] &&
            context->values[ATTENUATION] < 0) {
            update_filter_weights(context);
        }
    }

#ifndef NDEBUG
    struct timespec t_0;

    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &t_0);
#endif

    /* Filter the input samples, using overlap-save fast
     * convolution. */

    /* Given an impulse response length of M and an FFT length of N,
     * each round of overlap-save calculates N - M + 1 "useful"
     * filtered samples.  If the size is larger than that, we leave
     * the rest for the next transfer call. */

    if (size > N - M + 1) {
        size = N - M + 1;
    }

    /* Shift input fragments into the buffer, right to left,
     * interleaving them in the order that FFTW expects. */

    memmove(context->input, context->input + size * C,
            (N - size) * C * sizeof(context->input[0]));

    for (j = 0 ; j < C ; j += 1) {
        const snd_pcm_channel_area_t *a = &src_areas[j];
        unsigned char *s;
        float *t;

        s = (unsigned char *)(a->addr) + (a->first + a->step * src_offset) / 8;
        t = context->input + (N - size) * C + j;

        for (i = 0 ; i < size ; i += 1) {
            t[i * C] = (*(float *)(s + i * a->step / 8));
        }
    }

    /* Filter the samples normally if we're compensating, otherwise
     * copy them straight to the output.  Note that a copy is
     * necessary, because we want to keep the input buffer filled, in
     * order to be able to switch between compensation and straight
     * attenuation on the fly. */

    if (context->values[COMPENSATE] &&
        context->values[ATTENUATION] < 0) {
        fftwf_execute(context->input_to_bins);

        for (i = 0 ; i < N / 2 + 1 ; i += 1) {
            for (j = 0 ; j < C ; j += 1) {
                context->bins[C * i + j] *= context->weights[i];
            }
        }

        fftwf_execute(context->bins_to_output);

        /* FFTW calculates an unnormalized FFT transform, so the filter
         * output, needs to be scaled by 1 / N, when we're compensating. */

        c = 1.0 / N;
    } else {
        j = (N - size) * C;

        memcpy(context->output + j,
               context->input + j,
               size * C * sizeof(context->input[0]));

        /* When not compensating, we scale the input samples, to
         * achieve straight, uncompensated attenuation of the
         * specified amount. */

        c = expf(0.05 * logf(10) * (float)context->values[ATTENUATION]);
    }

    /* Copy the filtered samples into ALSA's buffer, de/interleaving
     * them as required. */

    for (j = 0 ; j < C ; j += 1) {
        const snd_pcm_channel_area_t *a = &dst_areas[j];
        float *s;
        unsigned char *t;

        s = context->output + (N - size) * C + j;
        t = (unsigned char *)(a->addr) + (a->first + a->step * dst_offset) / 8;

        for (i = 0 ; i < size ; i += 1) {
            *(float *)(t + i * a->step / 8) = s[i * C] * c;
        }
    }

#ifndef NDEBUG
    {
        struct timespec t;

        clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &t);
        context->benchmarks[0] += ((t.tv_sec - t_0.tv_sec) * 1000000000 +
                                   (t.tv_nsec - t_0.tv_nsec));
        context->benchmarks[1] += size;
    }
#endif

    return size;
}

static int init_callback(snd_pcm_extplug_t *ext)
{
#ifdef WITH_THREADS
    struct context *context = (struct context *)ext->private_data;

    if (!fftwf_init_threads()) {
        SNDERR("Could not initialize threading");
        return -1;
    }

    fftwf_plan_with_nthreads(context->threads);

    TRACE("Planning FFTs using %d threads.\n", context->threads);
#endif

    return 0;
}

static int hw_params_callback(snd_pcm_extplug_t *ext, snd_pcm_hw_params_t* params)
{
    struct context *context = (struct context *)ext->private_data;

    if(context->fft_length == 0) {
        /* Find the best FFT length based on maximum possible period size */

        snd_pcm_uframes_t psize;
        int ret, dir;
        if((ret = snd_pcm_hw_params_get_period_size_max(params, &psize, &dir)) < 0 || dir == 1) {
            SNDERR("could not query max period size");
            return -EINVAL;
        }

        /* Choose smallest power of two N such that N ≥ L+M-1 */
        context->fft_length = 1 << (int)ceilf(log2f(context->impulse_length + psize - 1));
        TRACE("Impulse length is %d and maximum period size is %ld, using optimal FFT size %d\n", context->impulse_length, psize, context->fft_length);
    }

    const unsigned int N = context->fft_length;
    unsigned int L;
    const int C = ext->channels;

#ifndef NDEBUG
    struct timespec t_0;

    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &t_0);
#endif

    /* Allocate buffers and create plans for the various FFT
     * transforms.  Filtering-related transforms operate on all
     * channels simultaneously and are set up to expect the samples
     * for each channel, interleaved in the same buffer, i.e. for 3
     * channels, the buffer would look like:

     |     1    |     2    |     3    |
     | c0 c1 c2 | c0 c1 c2 | c0 c1 c2 | ...

    */

    int plan_flags = context->wisdom_path ? FFTW_MEASURE : FFTW_ESTIMATE;
    context->interpolation_length = L = 1 << ((int)(ceilf(log2f((float)ext->rate / 2))));

    context->input = (float *)fftwf_malloc(sizeof(float) * N * C);
    context->output = (float *)fftwf_malloc(sizeof(float) * N * C);
    context->bins = (fftwf_complex *)fftwf_malloc(
        sizeof(fftwf_complex) * (N / 2 + 1) * C);
    context->weights = (fftwf_complex *)fftwf_malloc(
        sizeof(fftwf_complex) * (N / 2 + 1));

    context->input_to_bins = fftwf_plan_many_dft_r2c(
        1, (const int []) {N}, C,
        context->input, NULL, C, 1,
        context->bins, NULL, C, 1,
        plan_flags);

    context->bins_to_output = fftwf_plan_many_dft_c2r(
        1, (const int []) {N}, C,
        context->bins, NULL, C, 1,
        context->output, NULL, C, 1,
        plan_flags);

    /* This is necessary, in order to avoid filtering uninitialized
     * data during the first few transfer callbacks, before the input
     * buffer is filled. */

    memset(context->input, 0, sizeof(float) * N * C);

    /* Design-related transforms follow.  Here, no interleaving is
     * necessary. */

    context->kernel = (float *)fftwf_malloc(sizeof(float) * L);
    context->spline = (fftwf_complex *)fftwf_malloc(
        sizeof(fftwf_complex) * (L / 2 + 1));

    context->kernel_to_weights = fftwf_plan_dft_r2c_1d(N,
                                                       context->kernel,
                                                       context->weights,
                                                       plan_flags);

    context->spline_to_kernel = fftwf_plan_dft_c2r_1d(L,
                                                      context->spline,
                                                      context->kernel,
                                                      plan_flags);

    /* No more planning needs to be done, save FFTW wisdom if required */
    if(context->wisdom_path && fftwf_export_wisdom_to_filename(context->wisdom_path) != 1) {
        SNDERR("Failed saving wisdom to %s, continuing anyway", context->wisdom_path);
    }

#ifndef NDEBUG
    {
        struct timespec t;

        clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &t);
        TRACE("Initialized filter ("
              "FFT length %d, design length %d), in %.2f ms.\n",
              N, L,
              (t.tv_sec - t_0.tv_sec) * 1e3 +
              (t.tv_nsec - t_0.tv_nsec) * 1e-6);
    }
#endif

    if (!context->input || !context->output ||
        !context->bins || !context->weights ||
        !context->kernel || !context->spline ||
        !context->input_to_bins || !context->bins_to_output ||
        !context->spline_to_kernel || !context->kernel_to_weights) {
        /* The close_callback will be called to free any of these that
         * did get allocated. */

        return -1;
    } else if (context->values[COMPENSATE] &&
               context->values[ATTENUATION] < 0) {
        update_filter_weights(context);
    }

#ifdef PLOT_FILTER
    TRACE("Dumped filter data.  Will now force an exit.\n%s", "");

    return -1;
#endif

    return 0;
}

static int close_callback(snd_pcm_extplug_t *ext) {
    struct context *context = (struct context *)ext->private_data;

    /* Free all FFTW-related resources. */

    fftwf_free(context->input);
    fftwf_free(context->output);
    fftwf_free(context->bins);
    fftwf_free(context->weights);

    fftwf_destroy_plan(context->input_to_bins);
    fftwf_destroy_plan(context->bins_to_output);

    fftwf_free(context->kernel);
    fftwf_free(context->spline);

    fftwf_destroy_plan(context->spline_to_kernel);
    fftwf_destroy_plan(context->kernel_to_weights);

    /* Close the card's ctl interface. */

    if (context->ctl) {
        snd_ctl_close(context->ctl);
    }

#ifndef NDEBUG
    TRACE("Filtered %lu samples in %.3f s (at a rate of %.0f ns/sample).\n",
          context->benchmarks[1],
          context->benchmarks[0] * 1e-9f,
          context->benchmarks[0] / (float)context->benchmarks[1]);
#endif

    free(context->prefix);
    free(context);

    return 0;
}

static snd_pcm_extplug_callback_t callbacks = {
    .transfer = transfer_callback,
    .init = init_callback,
    .hw_params = hw_params_callback,
    .close = close_callback,
};

static int configure_card_ctl(struct context *context, int uninstall)
{
    snd_ctl_t *ctl;

    snd_pcm_info_t *pcm_info;
    snd_ctl_elem_id_t *id;
    snd_ctl_elem_info_t *ctl_info;

    char ctl_name[16];
    int e, card, device, subdevice;
    unsigned int i;

    /* Try to find the hardware card this PCM is connected to. */

    snd_pcm_info_alloca(&pcm_info);

    if (snd_pcm_info(context->ext.pcm, pcm_info) < 0) {
        return -1;
    }

    if ((card = snd_pcm_info_get_card(pcm_info)) < 0) {
        SNDERR("Could not get a card number from PCM");
        return -1;
    }

    /* Get a couple of parameters. */

    device = snd_pcm_info_get_device(pcm_info);
    subdevice = snd_pcm_info_get_subdevice(pcm_info);

    /* Open the ctl. */

    snprintf(ctl_name, sizeof(ctl_name), "hw:%d", card);

    if (snd_ctl_open(&ctl, ctl_name, SND_CTL_NONBLOCK) < 0) {
        return -1;
    }

    /* Iterate through all our controls and try to install them into
     * (or uninstall them from) the card's mixer. */

    snd_ctl_elem_id_alloca(&id);
    snd_ctl_elem_info_alloca(&ctl_info);

    for (i = 0 ; i < CONTROLS_COUNT ; i += 1) {
        snd_ctl_elem_type_t type;

        /* Figure out the element type. */

        if (CONTROLS[i].range[0] == 0 &&
            CONTROLS[i].range[1] == 1) {
            type = SND_CTL_ELEM_TYPE_BOOLEAN;
        } else {
            type = SND_CTL_ELEM_TYPE_INTEGER;
        }

        /* Configure the control element. */

        snd_ctl_elem_id_set_interface(id, SND_CTL_ELEM_IFACE_MIXER);
        if (!context->prefix) {
            snd_ctl_elem_id_set_name(id, CONTROLS[i].name);
        } else {
            char s[strlen(context->prefix) + strlen(CONTROLS[i].name) + 1];

            strcpy(s, context->prefix);
            strcat(s, CONTROLS[i].name);

            snd_ctl_elem_id_set_name(id, s);
        }

        snd_ctl_elem_id_set_index(id, 0);
        snd_ctl_elem_id_set_device(id, device);
        snd_ctl_elem_id_set_subdevice(id, subdevice);

        /* See if anything similar like this already exists. */

        snd_ctl_elem_info_set_id(ctl_info, id);
        if ((e = snd_ctl_elem_info(ctl, ctl_info)) < 0) {
            if (e != -ENOENT) {
                SNDERR("Cannot get info for ctl '%s'", ctl_name);
                snd_ctl_close(ctl);
                return -1;
            }

            /* The control doesn't exist; if uninstalling, we're done,
             * otherwise continue to install the control below. */

            if (uninstall) {
                continue;
            }
        } else {
            /* The control, or at least some control by the same name,
             * already exist. */

            if (snd_ctl_elem_info_is_user(ctl_info) == 0) {
                /* The existing control is a hardware control. Nothing
                 * much we can do, at this point. */

                SNDERR("Hardware control exists");
                snd_ctl_close(ctl);
                return -1;
            } else {
                /* The existing control is a user control, probably
                 * created by us, during a previous run.  See if its
                 * configuration matches our needs and reuse it,
                 * unless we're uninstalling. */

                if (!uninstall &&
                    snd_ctl_elem_info_get_count(ctl_info) == 1 &&
                    snd_ctl_elem_info_get_type(ctl_info) == type &&
                    ((type == SND_CTL_ELEM_TYPE_BOOLEAN) ||
                     ((snd_ctl_elem_info_get_min(ctl_info) ==
                       CONTROLS[i].range[0]) &&
                      (snd_ctl_elem_info_get_max(ctl_info) ==
                       CONTROLS[i].range[1])))) {

                    snd_ctl_elem_value_t *v;

                    /* Looks like we'll be using the control.  Read
                     * its value and update the filter parameters. */

                    snd_ctl_elem_value_alloca(&v);
                    snd_ctl_elem_value_set_id(v, id);

                    if (snd_ctl_elem_read(ctl, v) < 0) {
                        snd_ctl_close(ctl);
                        return -1;
                    }

                    context->values[i] = snd_ctl_elem_value_get_integer(v, 0);

                    TRACE("Using existing user control with name '%s'.\n",
                          snd_ctl_elem_id_get_name(id));

                    /* We're done with this control; move on. */

                    continue;
                }

                /* The control seems to have the wrong configuration;
                 * remove it. */

                if (snd_ctl_elem_remove(ctl, id) < 0) {
                    SNDERR("Could not remove existing, unsuitable "
                           "user control");
                    snd_ctl_close(ctl);
                    return -1;
                }

                TRACE("Removed existing user control with name '%s'.\n",
                      snd_ctl_elem_id_get_name(id));

                /* If uninstalling, we're done, otherwise continue to
                 * reinstall the control below. */

                if (uninstall) {
                    continue;
                }
            }
        }

        /* Inject a control, into the card's mixer. */

        if (type == SND_CTL_ELEM_TYPE_BOOLEAN) {
            if (snd_ctl_add_boolean_elem_set(ctl, ctl_info, 1, 1) < 0) {
                snd_ctl_close(ctl);
                return -1;
            }
        } else {
            SNDRV_CTL_TLVD_DECLARE_DB_MINMAX(tlv,
                                             CONTROLS[i].range[0] * 100,
                                             CONTROLS[i].range[1] * 100);


            if (snd_ctl_add_integer_elem_set(ctl, ctl_info, 1, 1,
                                             CONTROLS[i].range[0],
                                             CONTROLS[i].range[1],
                                             1) < 0 ||
                snd_ctl_elem_tlv_write(ctl, id, tlv) < 0) {
                snd_ctl_close(ctl);
                return -1;
            }
        }

        TRACE("Added user control with name '%s'.\n",
              snd_ctl_elem_id_get_name(id));

        /* Configure the control's dB range metadata and initialize
         * its value. */

        {
            snd_ctl_elem_value_t *v;
            snd_ctl_elem_value_alloca(&v);

            snd_ctl_elem_value_set_id(v, id);
            snd_ctl_elem_value_set_integer(v, 0, context->values[i]);

            if (snd_ctl_elem_write(ctl, v) < 0 ||
                snd_ctl_elem_unlock(ctl, id) < 0) {
                snd_ctl_close(ctl);
                return -1;
            }
        }
    }

    /* If we're uninstalling, we're done; close the mixer. */

    if (uninstall) {
        snd_ctl_close(ctl);
        return 0;
    } else {
        /* Install a handler, to monitor subsequent changes to the
         * control's settings. */

        if ((e = snd_ctl_subscribe_events(ctl, 1)) < 0) {
            SNDERR("Could not subscribe for ctl events");
            snd_ctl_close(ctl);
            context->ctl = NULL;
            return -1;
        } else {
            context->ctl = ctl;
            return 0;
        }
    }
}

SND_PCM_PLUGIN_DEFINE_FUNC(loudness)
{
    snd_config_iterator_t i, next;
    struct context *context;
    snd_config_t *slave = NULL;
    long int reference = 82, attenuation = -10, compensate = 1;
    long int impulse_length = 4096, fft_length = 0;

#ifdef WITH_THREADS
    long int threads = 1;
#endif

    int install = 0, uninstall = 0;
    const char *prefix = NULL, *wisdom_path = NULL;

    int e;

    /* Parse plugin configuration. */

    snd_config_for_each(i, next, conf) {
        snd_config_t *n = snd_config_iterator_entry(i);
        const char *id;

        if (snd_config_get_id(n, &id) < 0 ||
            strcmp(id, "comment") == 0 ||
            strcmp(id, "type") == 0 ||
            strcmp(id, "hint") == 0) {
            continue;
        }

        if (strcmp(id, "slave") == 0) {
            slave = n;
            continue;
        }

        if (strcmp(id, "controls") == 0) {
            snd_config_iterator_t i, next;
            snd_config_for_each(i, next, n) {
                snd_config_t *n = snd_config_iterator_entry(i);
                const char *id;

                if (snd_config_get_id(n, &id) < 0) {
                    continue;
                }

                if (strcmp(id, "install") == 0) {
                    install = snd_config_get_bool(n);
                    continue;
                }

                if (strcmp(id, "uninstall") == 0) {
                    uninstall = snd_config_get_bool(n);
                    continue;
                }

                if (strcmp(id, "prefix") == 0) {
                    snd_config_get_string(n, &prefix);

                    continue;
                }

                SNDERR("Unknown field controls.%s", id);
                return -EINVAL;
            }

            continue;
        }

        if (strcmp(id, "reference") == 0) {
            snd_config_get_integer(n, &reference);

            if (reference < CONTROLS[REFERENCE].range[0] ||
                reference > CONTROLS[REFERENCE].range[1]) {
                SNDERR("Reference level must be between %ddB and %ddB",
                       CONTROLS[REFERENCE].range[0],
                       CONTROLS[REFERENCE].range[1]);
                return -EINVAL;
            }

            continue;
        }

        if (strcmp(id, "attenuation") == 0) {
            snd_config_get_integer(n, &attenuation);

            if (attenuation < CONTROLS[ATTENUATION].range[0] ||
                attenuation > CONTROLS[ATTENUATION].range[1]) {
                SNDERR("Attenuation level must be between %ddB and %ddB",
                       CONTROLS[ATTENUATION].range[0],
                       CONTROLS[ATTENUATION].range[1]);
                return -EINVAL;
            }

            continue;
        }

        if (strcmp(id, "window") == 0) {
            snd_config_get_integer(n, &impulse_length);

            if (impulse_length < 1024) {
                SNDERR("Window length must not be lower than 1024");
                return -EINVAL;
            }

            continue;
        }

        if (strcmp(id, "fft") == 0) {
            snd_config_get_integer(n, &fft_length);

            if (fft_length != 0 && fft_length < 2 * impulse_length) {
                SNDERR("FFT length must be at least double the impulse length");
                return -EINVAL;
            }

            continue;
        }

        if (strcmp(id, "threads") == 0) {
#ifdef WITH_THREADS
            snd_config_get_integer(n, &threads);

            if (threads < 1) {
                SNDERR("Cannot use fewer than one threads");
                return -EINVAL;
            }
#endif

            continue;
        }

        if(strcmp(id, "wisdom_path") == 0) {
            const char* s;
            snd_config_get_string(n, &s);
            wisdom_path = strdup(s);
            if(wisdom_path == NULL) {
                SYSERR("strdup");
                return -EINVAL;
            }
            continue;
        }

        SNDERR("Unknown field %s", id);
        return -EINVAL;
    }

    if (install && uninstall) {
        SNDERR("Cannot both install and uninstall controls");
        return -EINVAL;
    }

    if (!slave) {
        SNDERR("No slave defined for loudness plugin");
        return -EINVAL;
    }

    /* All seems to be in order; allocate and initialize a new
     * filtering context. */

    if ((context = calloc(1, sizeof(*context))) == NULL) {
        return -ENOMEM;
    }

    /* Initialize the filtering parameters, with the values from the
     * configuration for now.  We'll update them from the mixer later,
     * if possible. */

#ifdef WITH_THREADS
    context->threads = threads;
#endif

    context->impulse_length = impulse_length;
    context->fft_length = fft_length;

    if (prefix) {
        context->prefix = strdup(prefix);
        if(context->prefix == NULL) {
            SYSERR("strdup");
            e = -EINVAL;
            goto cleanup;
        }
    } else {
        context->prefix = NULL;
    }

    context->values[REFERENCE] = reference;
    context->values[ATTENUATION] = attenuation;
    context->values[COMPENSATE] = compensate;

    context->ext.version = SND_PCM_EXTPLUG_VERSION;
    context->ext.name = "Loudness compensation";
    context->ext.callback = &callbacks;
    context->ext.private_data = context;

    /* Try loading FFTW wisdom. */

    context->wisdom_path = wisdom_path;
    if(context->wisdom_path && fftwf_import_wisdom_from_filename(context->wisdom_path) != 1) {
        SNDERR("Failed loading wisdom from %s, continuing anyway", context->wisdom_path);
    }

    /* Create the plugin. */

    if ((e = snd_pcm_extplug_create(&context->ext, name, root,
                                    slave, stream, mode)) < 0) {
        goto cleanup;
    }

    /* Set PCM constraints. */

    if ((e = snd_pcm_extplug_set_param(&context->ext,
                                       SND_PCM_EXTPLUG_HW_FORMAT,
                                       SND_PCM_FORMAT_FLOAT)) < 0) {
        goto cleanup;
    }

    if ((e = snd_pcm_extplug_set_slave_param(&context->ext,
                                             SND_PCM_EXTPLUG_HW_FORMAT,
                                             SND_PCM_FORMAT_FLOAT)) < 0) {
        goto cleanup;
    }

    /* Try to open the card's ctl and install (or uninstall, if so
     * requested) our controls. */

    if ((install || uninstall) &&
        configure_card_ctl(context, uninstall) < 0) {
        SNDERR("Could not (un)install controls");
    }

    *pcmp = context->ext.pcm;

    return 0;

cleanup:
    if(context->prefix) free((void*)context->prefix);
    if(context->wisdom_path) free((void*)context->wisdom_path);
    free(context);
    return e;
}

SND_PCM_PLUGIN_SYMBOL(loudness)
