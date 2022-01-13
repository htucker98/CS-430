public class Pixel {
    public int color;


    public Pixel() {
        this.color = 0;
    }

    public int getColor() {
        return this.color;
    }

    public void setBlack() {
        this.color = 1;
    }

    public boolean isBlack() {
        if (this.color == 1) {
            return true;
        } else {
            return false;
        }
    }
}

