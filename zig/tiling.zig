const std = @import("std");

fn tanh(x: f32) f32 {
    const exp_x = std.math.exp(x);
    const exp_minus_x = std.math.exp(-x);
    return (exp_x - exp_minus_x) / (exp_x + exp_minus_x);
}

pub fn rampedSin(omega: f32, width: f32, delay: f32, q: comptime_int) [q]f32 {
    var result = std.mem.zeroes([q]f32);
    var t: f32 = undefined;
    for (0..q) |i| {
        t = @as(f32, @floatFromInt(i)) * omega;
        result[i] = ((1 + tanh(t / width - delay)) / 2.0) * std.math.sin(t);
    }
    return result;
}

pub fn newFields(k: usize, m: usize, n: usize, p: usize, c: usize) ![]f32 {
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    const allocator = gpa.allocator();
    const K = k * m * n * p * c;
    var EH: []f32 = try allocator.alloc(f32, K);
    for (0..K) |q| {
        EH[q] = 0.0;
    }
    return EH;
}

pub fn writeBytes(arr: []const u8, filename: []const u8) !void {
    const file = try std.fs.cwd().createFile(filename, .{ .read = true });
    try file.writeAll(arr);
    file.close();
}

pub fn main() !void {
    const base: usize = 8;
    const height: usize = 4;
    const M: usize = 160;
    const N: usize = 96;
    const P: usize = 64;
    const Q: usize = 96;
    const Sc: f32 = 0.99 / std.math.sqrt(3.0);
    const Sx: usize = M / 2;
    const Sy: usize = N / 2;
    const Sz: usize = P / 2;
    const St: [Q]f32 = rampedSin(0.3, 5.0, 3.0, Q);
    const maxOffsetIdx = (M / base) * (N / base) * (P / base);

    var c1: usize = undefined;
    var c2: usize = undefined;
    var i: usize = undefined;
    var j0: isize = undefined;
    var j1: isize = undefined;
    var j2: isize = undefined;
    var m1: usize = undefined;
    var m2: usize = undefined;
    var m: usize = undefined;
    var mm: usize = undefined;
    var n0: usize = undefined;
    var n2: usize = undefined;
    var n: usize = undefined;
    var nn: usize = undefined;
    var p0: usize = undefined;
    var p1: usize = undefined;
    var p: usize = undefined;
    var pp: usize = undefined;
    var q0: usize = undefined;
    var q1: usize = undefined;
    var s: isize = undefined;
    var xmax: usize = undefined;
    var xmin: usize = undefined;
    var xoffset: usize = undefined;
    var ymax: isize = undefined;
    var ymin: isize = undefined;
    var yoffset: usize = undefined;
    var zmax: isize = undefined;
    var zmin: isize = undefined;
    var zoffset: usize = undefined;
    var _n: isize = undefined;
    var _p: isize = undefined;
    var xoffsetidx: usize = undefined;
    var yoffsetidx: usize = undefined;
    var zoffsetidx: usize = undefined;

    const EH = try newFields(2, M, N, P, 3);
    //EH[(height - 2) * M * N * P * 3 + Sx * N * P * 3 + Sy * P * 3 + Sz * 3 + 2] = 1.0; // E-excitation
    const startTime = std.time.milliTimestamp();
    for (0..2 * Q / height) |_i| {
        i = _i * height;
        for (0..2) |mvx| {
            for (0..maxOffsetIdx) |offsetidx| {
                xoffsetidx = offsetidx / (N / base) / (P / base);
                yoffsetidx = offsetidx / (P / base);
                zoffsetidx = @mod(offsetidx, (P / base));
                xoffset = xoffsetidx * base;
                yoffset = yoffsetidx * base;
                zoffset = zoffsetidx * base;
                for (0..height) |h| {
                    s = 1 - @as(isize, @intCast(2 * (h % 2)));
                    q0 = (h % 2) * M * N * P * 3;
                    q1 = ((2 + h - 1) % 2) * M * N * P * 3;
                    xmin = xoffset + (if (mvx == 0) h / 2 else base - 1 - (h + 1) / 2);
                    xmax = xoffset + (if (mvx == 0) base - 1 - (h + 1) / 2 else base + h / 2);
                    ymin = @as(isize, @intCast(yoffset)) - @as(isize, @intCast((h + 1) / 2));
                    ymax = @as(isize, @intCast(yoffset + base - (h + 1) / 2));
                    zmin = @as(isize, @intCast(zoffset)) - @as(isize, @intCast((h + 1) / 2));
                    zmax = @as(isize, @intCast(zoffset + base - (h + 1) / 2));
                    for (xmin..xmax) |_m| {
                        m = @mod(_m, M);
                        _n = ymin;
                        while (_n < ymax) : (_n += 1) {
                            n = @as(usize, @intCast(@mod(_n, @as(isize, @intCast(N)))));
                            _p = zmin;
                            while (_p < zmax) : (_p += 1) {
                                p = @as(usize, @intCast(@mod(_p, @as(isize, @intCast(P)))));
                                for (0..3) |c0| {
                                    c1 = @mod(c0 + 1, 3);
                                    c2 = @mod(c0 + 2, 3);
                                    j0 = s * @as(isize, @intFromBool(c0 == 0));
                                    j1 = s * @as(isize, @intFromBool(c0 == 1));
                                    j2 = s * @as(isize, @intFromBool(c0 == 2));
                                    mm = m * N * P * 3;
                                    nn = n * P * 3;
                                    pp = p * 3;
                                    m2 = ((m + @as(usize, @intCast(@as(isize, @intCast(M)) - j2))) % M) * N * P * 3;
                                    n0 = ((n + @as(usize, @intCast(@as(isize, @intCast(N)) - j0))) % N) * P * 3;
                                    p1 = ((p + @as(usize, @intCast(@as(isize, @intCast(P)) - j1))) % P) * 3;
                                    m1 = ((m + @as(usize, @intCast(@as(isize, @intCast(M)) - j1))) % M) * N * P * 3;
                                    n2 = ((n + @as(usize, @intCast(@as(isize, @intCast(N)) - j2))) % N) * P * 3;
                                    p0 = ((p + @as(usize, @intCast(@as(isize, @intCast(P)) - j0))) % P) * 3;
                                    EH[q0 + mm + nn + pp + c0] += Sc * (EH[q1 + mm + nn + pp + c2] - EH[q1 + m2 + n0 + p1 + c2] - EH[q1 + mm + nn + pp + c1] + EH[q1 + m1 + n2 + p0 + c1]);
                                }
                                if ((h % 2 == 0) and (m == Sx) and (n == Sy) and (p == Sz)) {
                                    EH[q0 + Sx * N * P * 3 + Sy * P * 3 + Sz * 3 + 2] += St[(i + h) / 2];
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    const endTime = std.time.milliTimestamp() - startTime;
    std.debug.print("execution time: {}ms\n\n", .{endTime});
    const EHbytes: [*]u8 = @ptrCast(EH);
    try writeBytes(EHbytes[0 .. 4 * EH.len], "output.bin");
}
