public class Line {
    public Point q;
    public Point r;


    public Line(Point q, Point r){
        this.q=q;
        this.r=r;
    }

    public Point getQ(){
        return q;
    }

    public Point getR(){
        return r;
    }

    public void setQ(Point q){
        this.q=q;
    }

    public void setR(Point r){
        this.r=r;
    }

}