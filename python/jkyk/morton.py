def morton2D(x: int, y: int) -> int:
    return (_part1by1(x) << 1) | _part1by1(y)


def imorton2D(index: int) -> tuple[int, int]:
    x = _compact1by1(index >> 1)
    y = _compact1by1(index)
    return x, y


def morton3D(x: int, y: int, z: int) -> int:
    return (_part1by2(x) << 2) | (_part1by2(y) << 1) | (_part1by2(z) << 0)


def imorton3D(index: int) -> tuple[int, int, int]:
    z = _compact1by2(index >> 0)
    y = _compact1by2(index >> 1)
    x = _compact1by2(index >> 2)
    return x, y, z


def _part1by1(n):
    n &= 0x0000FFFF
    n = (n | (n << 8)) & 0x00FF00FF
    n = (n | (n << 4)) & 0x0F0F0F0F
    n = (n | (n << 2)) & 0x33333333
    n = (n | (n << 1)) & 0x55555555
    return n


def _compact1by1(n):
    n &= 0x55555555
    n = (n ^ (n >> 1)) & 0x33333333
    n = (n ^ (n >> 2)) & 0x0F0F0F0F
    n = (n ^ (n >> 4)) & 0x00FF00FF
    n = (n ^ (n >> 8)) & 0x0000FFFF
    return n


def _part1by2(n):
    n &= 0x000003FF  # only 10 bits
    n = (n | (n << 16)) & 0x030000FF
    n = (n | (n << 8)) & 0x0300F00F
    n = (n | (n << 4)) & 0x030C30C3
    n = (n | (n << 2)) & 0x09249249
    return n


def _compact1by2(n):
    n &= 0x09249249
    n = (n ^ (n >> 2)) & 0x030C30C3
    n = (n ^ (n >> 4)) & 0x0300F00F
    n = (n ^ (n >> 8)) & 0x030000FF
    n = (n ^ (n >> 16)) & 0x000003FF
    return n
