package com.langqi.fakegps;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;

/**
 * 根据绕操场跑道，按目标总里程生成随机多圈路线。
 * 算法移植自独立工具 KmlTrackGenerator；文件读写与起点计数持久化由调用方负责。
 */
public class TrackGenerator {

    /** 地球平均半径（米） */
    static final double R = 6371000.0;

    /** 要求的里程余量（米） */
    public static final double MARGIN_M = 100.0;

    /** 圈数检测垂线段的半长（米），总长 20 米 */
    static final double CROSS_LINE_HALF_M = 10.0;

    /** 每纬度对应的大致距离（米） */
    static final double M_PER_DEG_LAT = 111320.0;

    /** 每一圈内部点的随机扰动标准差（米），起点除外 */
    static final double JITTER_SIGMA_M = 2.0;

    /** 起点偏移每个区间的宽度（米） */
    static final double START_OFFSET_STEP_M = 10.0;

    /** 起点偏移区间的个数（覆盖 0 ~ 50 米） */
    public static final int START_OFFSET_BUCKETS = 5;

    /** 一个经纬度点（顺序与 KML 坐标一致：lon, lat） */
    public static class Point {
        public final double lon;
        public final double lat;

        public Point(double lon, double lat) {
            this.lon = lon;
            this.lat = lat;
        }
    }

    /**
     * 根据参考跑道生成随机多圈路线。
     *
     * @param trackPoints  参考跑道坐标（应为闭合多圈）
     * @param totalMeters  目标总里程（米，不含余量）
     * @param marginMeters 余量（米）
     * @param startBucket  起点偏移区间序号 0..START_OFFSET_BUCKETS-1
     * @param rng          随机源
     * @return 生成的路径；若参考跑道无法识别为闭合多圈则返回 null
     */
    public static List<Point> generate(List<Point> trackPoints, double totalMeters,
                                       double marginMeters, int startBucket, Random rng) {
        if (trackPoints == null || trackPoints.size() < 2) {
            return null;
        }
        double totalLength = polylineLength(trackPoints);
        List<Integer> crossings = findCrossings(trackPoints);
        int laps = crossings.size();
        if (laps < 1) {
            return null;
        }
        double singleLapPerimeter = totalLength / laps;
        if (singleLapPerimeter <= 0) {
            return null;
        }

        List<Point> singleLap = extractSingleLap(trackPoints, crossings);
        rescaleToPerimeter(singleLap, singleLapPerimeter);

        double startOffset = startBucket * START_OFFSET_STEP_M + rng.nextDouble() * START_OFFSET_STEP_M;
        singleLap = cutLoopAt(singleLap, startOffset);

        double target = totalMeters + marginMeters;
        int fullLaps = (int) (target / singleLapPerimeter);
        double remainder = target - fullLaps * singleLapPerimeter;

        return buildPath(singleLap, singleLapPerimeter, fullLaps, remainder, rng);
    }

    /** Haversine 距离（米） */
    static double haversine(Point a, Point b) {
        double lat1 = Math.toRadians(a.lat), lat2 = Math.toRadians(b.lat);
        double dLat = Math.toRadians(b.lat - a.lat);
        double dLon = Math.toRadians(b.lon - a.lon);
        double h = Math.sin(dLat / 2) * Math.sin(dLat / 2)
                + Math.cos(lat1) * Math.cos(lat2) * Math.sin(dLon / 2) * Math.sin(dLon / 2);
        return 2 * R * Math.asin(Math.sqrt(h));
    }

    /** 折线总长度（按相邻点累加） */
    static double polylineLength(List<Point> pts) {
        double sum = 0;
        for (int i = 0; i < pts.size() - 1; i++) {
            sum += haversine(pts.get(i), pts.get(i + 1));
        }
        return sum;
    }

    /** 闭环长度（假定首尾相连） */
    static double closedLoopLength(List<Point> loop) {
        double sum = 0;
        for (int i = 0; i < loop.size() - 1; i++) {
            sum += haversine(loop.get(i), loop.get(i + 1));
        }
        return sum;
    }

    /** 检测轨迹穿越起点垂线段的段索引列表（穿越一次计一圈） */
    static List<Integer> findCrossings(List<Point> pts) {
        List<Integer> crossings = new ArrayList<>();
        if (pts.size() < 2) {
            return crossings;
        }
        Point p0 = pts.get(0);
        Point p1 = pts.get(1);
        double lon0 = p0.lon, lat0 = p0.lat;

        double[] P1 = toLocal(p1, lon0, lat0);
        double mx = P1[0] / 2, my = P1[1] / 2;
        double nx = -P1[1], ny = P1[0];
        double nl = Math.hypot(nx, ny);
        nx /= nl;
        ny /= nl;
        double[] A = {mx - nx * CROSS_LINE_HALF_M, my - ny * CROSS_LINE_HALF_M};
        double[] B = {mx + nx * CROSS_LINE_HALF_M, my + ny * CROSS_LINE_HALF_M};

        for (int i = 0; i < pts.size() - 1; i++) {
            double[] q1 = toLocal(pts.get(i), lon0, lat0);
            double[] q2 = toLocal(pts.get(i + 1), lon0, lat0);
            if (segmentsIntersect(q1, q2, A, B)) {
                crossings.add(i);
            }
        }
        return crossings;
    }

    /** 线段相交判断（含端点共线情况） */
    static boolean segmentsIntersect(double[] p1, double[] p2, double[] p3, double[] p4) {
        double d1 = cross(p3, p4, p1);
        double d2 = cross(p3, p4, p2);
        double d3 = cross(p1, p2, p3);
        double d4 = cross(p1, p2, p4);
        if (((d1 > 0 && d2 < 0) || (d1 < 0 && d2 > 0)) &&
                ((d3 > 0 && d4 < 0) || (d3 < 0 && d4 > 0))) {
            return true;
        }
        if (d1 == 0 && onSegment(p3, p1, p4)) return true;
        if (d2 == 0 && onSegment(p3, p2, p4)) return true;
        if (d3 == 0 && onSegment(p1, p3, p2)) return true;
        if (d4 == 0 && onSegment(p1, p4, p2)) return true;
        return false;
    }

    static double cross(double[] o, double[] a, double[] b) {
        return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0]);
    }

    static boolean onSegment(double[] p, double[] q, double[] r) {
        return q[0] <= Math.max(p[0], r[0]) && q[0] >= Math.min(p[0], r[0]) &&
                q[1] <= Math.max(p[1], r[1]) && q[1] >= Math.min(p[1], r[1]);
    }

    /** 用检测到的穿越点提取第一圈作为单圈闭环 */
    static List<Point> extractSingleLap(List<Point> pts, List<Integer> crossings) {
        int a = crossings.get(0);
        int b = crossings.size() >= 2 ? crossings.get(1) : pts.size() - 1;
        List<Point> loop = new ArrayList<>();
        for (int i = a; i <= b; i++) {
            loop.add(pts.get(i));
        }
        loop.add(pts.get(a));
        return loop;
    }

    /** 以起点为锚点，等比缩放闭环到目标周长（保持起点不变） */
    static void rescaleToPerimeter(List<Point> loop, double targetPerimeter) {
        scaleFromAnchor(loop, targetPerimeter);
    }

    /** 把闭环从给定弧长偏移处切开，重新组成新的闭环（新的起点） */
    static List<Point> cutLoopAt(List<Point> loop, double arcOffset) {
        double perimeter = closedLoopLength(loop);
        double target = arcOffset % perimeter;
        if (target < 1e-9) {
            return new ArrayList<>(loop);
        }
        double acc = 0;
        int k = 0;
        for (int i = 0; i < loop.size() - 1; i++) {
            double d = haversine(loop.get(i), loop.get(i + 1));
            if (acc + d >= target) {
                k = i;
                break;
            }
            acc += d;
        }
        Point cut = interpolate(loop.get(k), loop.get(k + 1),
                (target - acc) / haversine(loop.get(k), loop.get(k + 1)));

        List<Point> out = new ArrayList<>();
        out.add(cut);
        for (int i = k + 1; i < loop.size() - 1; i++) {
            out.add(loop.get(i));
        }
        for (int i = 0; i <= k; i++) {
            out.add(loop.get(i));
        }
        out.add(cut);
        return out;
    }

    /** 生成一圈独立的扰动圈，并缩放到目标周长 */
    static List<Point> generateLap(List<Point> loop, double targetPerimeter, Random rng) {
        List<Point> noisy = jitterLoop(loop, rng, JITTER_SIGMA_M);
        scaleFromAnchor(noisy, targetPerimeter);
        return noisy;
    }

    /** 对一圈做逐点随机抖动，起点（首尾点）保持不变 */
    static List<Point> jitterLoop(List<Point> loop, Random rng, double sigma) {
        Point anchor = loop.get(0);
        double lon0 = anchor.lon, lat0 = anchor.lat;
        List<Point> out = new ArrayList<>(loop.size());
        out.add(anchor);
        for (int i = 1; i < loop.size() - 1; i++) {
            double[] p = toLocal(loop.get(i), lon0, lat0);
            double nx = p[0] + rng.nextGaussian() * sigma;
            double ny = p[1] + rng.nextGaussian() * sigma;
            out.add(fromLocal(nx, ny, lon0, lat0));
        }
        out.add(anchor);
        return out;
    }

    /** 以起点（首点）为锚点等比缩放，使闭环周长精确等于目标值，起点不动 */
    static void scaleFromAnchor(List<Point> loop, double targetPerimeter) {
        double current = closedLoopLength(loop);
        if (current <= 0) return;
        double scale = targetPerimeter / current;
        Point anchor = loop.get(0);
        double lon0 = anchor.lon, lat0 = anchor.lat;
        for (int i = 1; i < loop.size() - 1; i++) {
            double[] p = toLocal(loop.get(i), lon0, lat0);
            loop.set(i, fromLocal(p[0] * scale, p[1] * scale, lon0, lat0));
        }
        loop.set(loop.size() - 1, anchor);
    }

    /** 组装完整路径：逐圈独立生成（每圈不同），再补足剩余里程 */
    static List<Point> buildPath(List<Point> loop, double lapPerimeter,
                                 int fullLaps, double remainder, Random rng) {
        List<Point> path = new ArrayList<>();
        for (int l = 0; l < fullLaps; l++) {
            List<Point> lap = generateLap(loop, lapPerimeter, rng);
            for (int i = 0; i < lap.size() - 1; i++) {
                path.add(lap.get(i));
            }
        }
        if (remainder > 1e-6) {
            List<Point> lap = generateLap(loop, lapPerimeter, rng);
            path.add(lap.get(0));
            double acc = 0;
            for (int i = 1; i < lap.size(); i++) {
                double d = haversine(lap.get(i - 1), lap.get(i));
                if (acc + d >= remainder) {
                    double frac = (remainder - acc) / d;
                    path.add(interpolate(lap.get(i - 1), lap.get(i), frac));
                    break;
                }
                path.add(lap.get(i));
                acc += d;
            }
        }
        return path;
    }

    static Point interpolate(Point a, Point b, double t) {
        return new Point(a.lon + (b.lon - a.lon) * t, a.lat + (b.lat - a.lat) * t);
    }

    /** 本地平面坐标（米），用于等比缩放与扰动 */
    static double[] toLocal(Point c, double lon0, double lat0) {
        double x = (c.lon - lon0) * M_PER_DEG_LAT * Math.cos(Math.toRadians(lat0));
        double y = (c.lat - lat0) * M_PER_DEG_LAT;
        return new double[]{x, y};
    }

    static Point fromLocal(double x, double y, double lon0, double lat0) {
        double lon = lon0 + x / (M_PER_DEG_LAT * Math.cos(Math.toRadians(lat0)));
        double lat = lat0 + y / M_PER_DEG_LAT;
        return new Point(lon, lat);
    }
}
