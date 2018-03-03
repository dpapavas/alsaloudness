function c = contour(SPL)
    [c_iso, f_iso] = iso226(SPL);
    f = [0.1 f_iso 22e3 100e3];

    c = csape(log(f), [1.75 * c_iso(1), c_iso, c_iso(1), 1.15 * c_iso(1)],
              'variational', [0;0]);
end
