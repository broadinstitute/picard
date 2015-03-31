package picard.illumina.quality;

public class Point {
    private int x;
    private int y;

    public Point(final int x, final int y) {
        this.x=x;
        this.y=y;
    }

    public int getX() { return x; }
    public int getY() { return y; }

    public double distance(Point p) {
        final int dx = p.x - x;
        final int dy = p.y - y;
        return Math.sqrt(dx*dx + dy*dy);
    }
}
