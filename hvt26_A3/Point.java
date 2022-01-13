public class Point{
    public int x;
    public int y;
    public double bitcode;



    public Point(int x, int y) {
        this.x =x;
        this.y =y;
        this.bitcode= 0;
    }

    public int getX(){
        return this.x;
    }

    public int getY(){
        return this.y;
    }

    public void setX(int x){
        this.x = x;
    }
    public void setY(int y) {
        this.y = y;
    }

    public boolean equals(Point q) {
        return x == q.x && y == q.y;
    }


}
