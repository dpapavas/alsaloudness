`alsaloudness` is a loudness-compensated volume control for ALSA.  It is based on the ISO 226 equal loudness contours, it is fast (thanks to FFTW for the most part) and should be able to run in real time, presenting a very small load to the CPU, even for the embedded systems often used in media player setups.

## Why would I want to use it? ##

There are plenty of sources describing the concepts of [loudness](https://en.wikipedia.org/wiki/Loudness) and [loudness compensation](https://en.wikipedia.org/wiki/Loudness_compensation) in detail, but the gist of it is this: the human ear has different sensitivity across the frequency spectrum, meaning that two tones of the same sound pressure, but different pitch, will not sound equally loud.  What's more if you reduce the level by the same amount, the perceived change in loudness won't be the same either.

The upshot of all this, is that as you lower the volume through a typical volume knob, which ultimately controls the level of the reproduced sound, you don't only get the intended reduction in volume; you also get an effective change in equalization, as bass and treble content is attenuated more than mid-range content.  As the level is progressively reduced further, larger and larger parts of the spectrum fall below the hearing threshold and the overall effect, is that the material sounds "thin".

A loudness-compensated volume control tries to rectify this situation, by suitably boosting the bass and treble parts of the spectrum, depending on the setting of the volume knob.

## I'm sold.  So, how do I use alsaloudness? ##

This involves a couple of steps.

### Compilation and installation ###

Compilation and installation of the module, can be achieved through the usual `make && sudo make install`, which, by default, will place the plugin in `/usr/local/lib`.  Once the plugin is installed, the per-user or global ALSA configuration (typically residing in `~/.asoundrc` and `/etc/asound.conf` respectively), will need to be amended as shown in the provided example `asoundrc` file.

On 64-bit systems running 32-bit applications, you must also install the 32-bit version of the plugin by using `make clean && make LIB32=yes && sudo make install LIB32=yes`.

If everything works out, the `loudness` PCM should now be available.  This can be set as the default PCM, if loudness compensation is desired for all audio output, otherwise only specific applications, such as MPD for instance, can be directed to use it.

One way to test whether the setup is working, is through the use of `aplay` as follows:

```bash
# aplay -Dloudness /path/to/some.wav
```

### Using the compensated controls ###

After the plugin is used for the first time, it will install a couple of new controls into the mixer of the sound card it plays through.  The default names of these controls, are `Attenuation` and `Reference`, but if your mixer already has controls by these names, they can be altered by specifying a prefix, as shown in the example `asoundrc` file.

To use these controls, the following strategy can be employed:

1. Set `Attenuation` to 0dB, effectively disabling loudness compensation and set your card's usual volume control, or the volume knob on your external amplifier or speakers, if that's what you normally use, at a level where the reproduced material "sounds right".  (Usually, as you raise the volume, the sound progressively gets "fuller", until it reaches a point which "sounds right".  From there on, further increase in volume will make the bass sound boomy and generally lead to a reduction in sound quality.)

2. Once the "proper" level is set, you should be able to reduce the reproduction volume to the desired listening level, by setting the `Attenuation` control appropriately, while, hopefully, maintaining the overall equalization constant.

Since this "proper" level will, generally, depend on the material, adjustments are likely to be necessary.  One can increase the compensation, i.e. further boost bass and treble content, while keeping the volume constant, by simultaneously increasing the "real", that is uncompensated, volume, while decreasing the `Attenuation` setting by the same amount.  The reverse is of course also possible, and would lead to a reduction in compensation.  These operations can easily be automated through a bit of scripting.

### Theory of operation (or, the good) ###

The theory behind all this, which also provides an explanation for the `Reference` control is as follows:

Although the human ear hears differently across the spectrum, this is typically taken care of when audio is produced.  If bass or treble needs to be boosted, because the human ear is less sensitive to these bands, the proper equalization will be applied during mastering, or when mixing live audio, or even when arranging acoustic instruments on the stage, until the desired listening experience is achieved.

One can therefore assume that audio, which has been produced at a certain reference level, will sound "right" when reproduced at the same level.  If then, the change in perceived equalization when listening at a lower level were known, one could, in theory, compensate for it, by boosting the parts of the spectrum that have been attenuated more strongly by the volume reduction, as required.

This is to some extent possible and `alsaloudness` attempts to accomplish this, by making use of so-called equal loudness contours, i.e. curves representing the level that a tone needs to have at each frequency, in order to have a certain loudness.  Such contours have been derived through psycho-acoustic experiments for a wide range of loudness levels and published as an ISO standard (ISO 226:2003, in its latest incarnation at the time of writing).  Given two such equal loudness curves, one for the reference level and one for the listening level, the desired compensation equalization curve can be simply derived, by subtracting the former from the latter.

The resulting equalization curves, as implemented (and produced) by `alsaloudness` can be seen below for several attenuation levels (relative to  a reference level of 90dB):

![Loudness compensation curves](./filterplots.png)

So, to get ideal equalization, one would have to know the reference level, at which the material was produced, set the `Reference` control accordingly, then set `Attenuation` to 0dB and adjust the volume until the same reference level is achieved.  After that, lowering the reproduction volume through the `Attenuation` control, should result in ideally compensated attenuation and hence no loss in sound fidelity.  In practice, the value for the `Reference` control, doesn't seem to make too much difference, so the default can probably be used most of the time.

### Operation in practice (or, the bad and the ugly) ###

In practice, several aspects of the operation laid out above are implemented in a less than ideal fashion:  For one, the equal loudness contours, that are the basis for compensation, are specified only up to 12.5kHz and need to be estimated above that.  Additionally, even in the range for which data are available, the contours are derived for listening conditions which will generally differ from the conditions of reproduction.

Add to that, the fact that the reference level is generally different depending on material and that it is generally unknown and it becomes evident that, whatever compensation is performed, will generally be less than ideal.  And even if it were ideal, there are probably other aspects connected to the volume of reproduction, that affect the listening experience, so that listening at lower volumes will never sound as good.

One can still ask whether listening at lower volumes can at least be improved though and hopefully `alsaloudness` is a step in the right direction.

### Further configuration options ###

There are a few more aspects of the plugin's configuration, that can be changed through the `asoundrc` file, although the default values should be fine in most cases.

Compensation is done through FIR filtering, implemented via FFT convolution, where the filter coefficients, are recalculated each time one of the parameters (i.e. `Attenuation` or `Reference`) is changed.  The default filter length is 4096 taps, which should work well in most cases, but this value can be changed though the `window` configuration variable, if necessary.  One possible reason to do so, might be to lower the load on the CPU, but it might also cause insufficient boosting of the bass frequencies (which is also the case for the default length of 4096, as is visible in the plots above, albeit at frequencies well below the threshold of hearing).

To examine the actual filter response calculated by the plugin, it should be recompiled with `make PLOT=yes && sudo make install`.  When built in this way, attempting to play through the plugin, will cause it to dump the filter parameters into a set of files in the current directory (named `parameters`, `sampled`, `realized`), along with a `gnuplot` script (named `filter.gnu`), which can be used to plot the generated filter response through the command `gnuplot -p filter.gnu`.

By default, the plugin determines the optimal value for the FFT size (`fft` parameter), based on the maximum period size requested by the application. This behaviour can be overridden by setting the `fft` parameter to any non-zero value. Good values will probably fall in the range of four to eight times the impulse length and be powers of two, or, at least, products of small prime factors.

Finally, the plugin is built with threading by default, which can be disabled, by compiling with `make THREADED=no && sudo make install`.  The number of used threads can be set via the `threads` variable in `asoundrc` and the default setting of 1, effectively disables threading.  The optimal number generally depends on the machine.  On desktop computers, using more than one threads, seems to be detrimental to performance, but using a thread count equal to the number of cores, seems to help on the Raspberry Pi 3.

To evaluate the filtering efficiency of some configuration, the plugin can be recompiled via `make BENCHMARK=yes && sudo make install`, which will cause it to print out timing information to the standard output.  Playing a small sound clip, via `aplay -d 10` for instance, should then allow evaluation of a configuration change.

### More drastic modifications ###

As already explained, the response of the compensation filter is based on equal loudness contours, derived from ISO 226 data.  Since this data is only available for part of the audible spectrum, some guesswork has been involved, in estimating the contours for higher frequencies.  If the resulting curves are deemed suspect, they can be changed relatively easily through a set of GNU Octave scripts, provided in the `/scripts` directory.

The contours themselves are interpolated from ISO 226 data points, inside the `contour.m` script, which is used by both `contours.m` and `validate.m`.  The former calculates and plots the interpolated contours and it also dumps a C header file, which can be compiled with the plugin.  The `validate.m` script, reads the files generated by the plugin, when built via `make PLOT=yes` and plots the realized filter response.  It also calculates the maximum deviation of certain parts of the filter computation from expected values, as a sanity check.

## License ##

`alsaloudness` is released under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.  See the file `LICENSE` for more details.
