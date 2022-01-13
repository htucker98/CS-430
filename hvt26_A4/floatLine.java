public class floatLine {

    public floatPoint q;
    public floatPoint r;


    public floatLine(floatPoint q, floatPoint r){
        this.q=q;
        this.r=r;
    }
    public floatPoint getQ(){
        return q;
    }

    public floatPoint getR(){
        return r;
    }

    public void setQ(floatPoint q){
        this.q=q;
    }

    public void setR(floatPoint r){
        this.r=r;
    }
}
