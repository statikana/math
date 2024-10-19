def max_magnitude(x, y, z):
    ax, ay, az = map(abs, (x, y, z))

    if ax > ay:
        # y x z
        # y z x
        if ax > az:
            return x
        else:
            return z
    else:
        if ay > az:
            return y  # y > x, y > z
        else:
            return z  # y > x, z > y
