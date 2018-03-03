pkg load splines

clf
clear

write = true;
SPL_ref = 90;
f_iso = nthargout (2, @iso226, 90);

pp_ref = contour(SPL_ref);

hold on

if write;
    fd = fopen("contours.h", "w");
    fprintf(fd, "#ifndef _CONTOURS_H\n#define _CONTOURS_H\n\nstatic struct {\n    float knots[33], coefficients[91][%d][%d];\n} contours = {\n    {\n        ",
            size(pp_ref.coefs));

    for j = 1:length(pp_ref.breaks);
        if mod(j, 6) == 0;
            fprintf(fd, "%.7f,\n        ", pp_ref.breaks(j));
        else
            fprintf(fd, "%.7f, ", pp_ref.breaks(j));
        end
    end

    fprintf(fd, "\n    },\n\n    {\n");
end

for SPL = 0 : 1 : 90;
    pp = contour(SPL);

    if write;
        fprintf(fd, "        {\n");

        for i = 1:size(pp_ref.coefs)(1);
            fprintf(fd, "            {");

            for j = 1:size(pp_ref.coefs)(2);
                fprintf(fd, "%.7f, ", pp.coefs(i, j));
            end

            fprintf(fd, "},\n");
        end

        fprintf(fd, "        },\n\n");
    end

    if mod(SPL_ref - SPL, 10) == 0;
        x = logspace(-1, 5, 1000);
        semilogx(x, ppval(pp, log(x)) - ppval(pp_ref, log(x)))
    end

    ylim ([-100 0])
    grid on;
end

if write;
    fprintf(fd, "    }\n};\n\n#endif");
    fclose(fd);
end
