import java.io.File;
import java.io.FileReader;
import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;



public class Image {
    public File ps;
    public File pbm;
    public String fileName;
    public ArrayList<Line> lines;
    public Pixel[][] pixels;
    public Window window;
    public int width;
    public int length;


    public Image(String[] args) throws FileNotFoundException{

        //default args
        String f = System.getProperty("user.dir") + "/" + "hw1.ps";
        int a = 0;
        int b = 0;
        int c = 499;
        int d = 499;
        String fileName = "hw1";

        // check if there is specific commandline arguments
        for(int i=0; i< args.length; i++){
            if(args[i].equals("-f")){
                f = System.getProperty("user.dir") + "/" + args[i+1];
                fileName = f.replace(".ps","");
            }
            if(args[i].equals("-a")){
                a = Integer.parseInt(args[i+1]);
            }
            if(args[i].equals("-b")){
                b = Integer.parseInt(args[i+1]);
            }
            if(args[i].equals("-c")){
                c = Integer.parseInt(args[i+1]);
            }
            if(args[i].equals("-d")){
                d = Integer.parseInt(args[i+1]);
            }
        }


        //check if ps file exists
        try {
            this.ps = new File(f);
            if (!this.ps.exists()) {
                throw new FileNotFoundException();
            }
        }
        catch (FileNotFoundException e) {
            e.printStackTrace();
            return;
        }

        //set fileName
        this.fileName = fileName;


        //intialize class attributes
        this.lines = new ArrayList<Line>();
        this.width = (c-a)+1;
        this.length = (d-b)+1;
        this.window = new Window(a, b, c, d);
        this.pixels = new Pixel[width][length];
        for (int i = 0; i <(this.width); i++) {
            for (int j = 0; j <(this.length); j++) {
                pixels[i][j] = new Pixel();
            }
        }
    
    }

    public void writePixel(int x, int y) {
        pixels[x][y].setBlack();
    }

    public void Scan() throws Exception {
        BufferedReader scanner = new BufferedReader(new FileReader(this.ps));
        //read file line by line
        String line = "";
        while (!line.contains("%%%END")) {
            line = scanner.readLine();
            //if line has the string "line" get the point values
            if (line.contains(" Line")) {
                String[] ints = line.split(" ");
                int x1 = Integer.parseInt(ints[0]);
                int y1 = Integer.parseInt(ints[1]);
                int x2 = Integer.parseInt(ints[2]);
                int y2 = Integer.parseInt(ints[3]);
                lines.add(new Line(new Point(x1, y1), new Point(x2, y2)));
            }
        }
    }

    public void bresenham() {
        // refrenced at; http://tech-algorithm.com/articles/drawing-line-using-bresenham-algorithm/
        for (Line line : this.lines) {
            int x = line.getQ().getX();
            int y = line.getQ().getY();
            int w = line.getR().getX() - line.getQ().getX();
            int h = line.getR().getY() - line.getQ().getY();
            int dx1 = 0;
            int dy1 = 0;
            int dx2 = 0;
            int dy2 = 0;
            if (w < 0) dx1 = -1;
            else if (w > 0) dx1 = 1;
            if (h < 0) dy1 = -1;
            else if (h > 0) dy1 = 1;
            if (w < 0) dx2 = -1;
            else if (w > 0) dx2 = 1;
            int longest = Math.abs(w);
            int shortest = Math.abs(h);
            if (!(longest > shortest)) {
                longest = Math.abs(h);
                shortest = Math.abs(w);
                if (h < 0) dy2 = -1;
                else if (h > 0) dy2 = 1;
                dx2 = 0;
            }
            int numerator = longest >> 1;
            for (int i = 0; i <= longest; i++) {
                this.writePixel(Math.abs(x-this.window.left), Math.abs(this.length-1)-(y-(this.window.bottom)));
                numerator += shortest;
                if (!(numerator < longest)) {
                    numerator -= longest;
                    x += dx1;
                    y += dy1;
                } else {
                    x += dx2;
                    y += dy2;
                }
            }
        }

    }

    public void assignBitCode(Point point) {

        if (point.getY() > this.window.top) {
            point.bitcode += 1000;
        }
        if (point.getY() < this.window.bottom) {
            point.bitcode += 100;
        }
        if (point.getX() > this.window.right) {
            point.bitcode += 10;
        }
        if (point.getX() < this.window.left) {
            point.bitcode += 1;
        }


    }

    public void cohenSutherland() {
        // keep index of lines in case it needs to be removed
        int i = 0;
        //list of outofbounds lines to remove
        ArrayList<Integer> outOfBounds = new ArrayList<>();

        for (Line line : this.lines) {
            double top = 1000;
            double bottom = 100;
            double right = 10;
            double left = 1;

            //assign bitcode to endpoints
            this.assignBitCode(line.getQ());
            this.assignBitCode(line.getR());

            // while loop to clip lines while any point on line is out of bounds
            while(true){

                //if both points in bounds
                if(line.getQ().bitcode == 0 && line.getR().bitcode==0){
                    break;
                }

                //if line is not visible in window
                if ((line.getQ().bitcode != 0 && line.getR().bitcode != 0) && ((line.getQ().getX() == line.getR().getX()) || (line.getQ().getY() == line.getR().getY()))) {
                    // don't write any pixels entire line is out of window
                    outOfBounds.add(i);
                    break;
                }

                else {
                    //if point Q is out of bounds
                    if (line.getQ().bitcode != 0) {
                        double codeOut;
                        int x = 0;
                        int y = 0;
                        int x1 = line.getQ().getX();
                        int x2 = line.getR().getX();
                        int y1 = line.getQ().getY();
                        int y2 = line.getR().getY();


                        codeOut = line.getQ().bitcode;

                        // Find intersection point

                        //If line intersects top of window
                        if (Math.round(codeOut / top) == 1) {
                            try {
                                x = x1 + (x2 - x1) * (window.top - y1) / (y2 - y1);
                            } catch (ArithmeticException arithmeticException){
                                x = x1;
                            }
                            y = window.top;
                        }
                        //If line intersects bottom of window
                        else if (Math.round(codeOut / bottom) == 1) {
                            try{
                                x = x1 + (x2 - x1) * (window.bottom - y1) / (y2 - y1);
                            } catch (ArithmeticException arithmeticException){
                                x = x1;
                            }
                            y = window.bottom;

                        }
                        //If line intersects right of window
                        else if (Math.round(codeOut / right) == 1) {
                            x = window.right;
                            try {
                                y = y1 + (y2 - y1) * (window.right - x1) / (x2 - x1);
                            } catch (ArithmeticException arithmeticException){
                                y =y1;
                            }
                        }
                        //If line intersects left of window
                        else if (codeOut / left == 1) {
                            x = window.left;
                            try {
                                y = y1 + (y2 - y1) * (window.left - x1) / (x2 - x1);
                            } catch (ArithmeticException arithmeticException){
                                y = y1;
                            }
                        }
                        Point p = new Point(x, y);
                        assignBitCode(p);
                        line.setQ(p);
                    }
                    if (line.getR().bitcode != 0) {
                        double codeOut;
                        int x = 0;
                        int y = 0;
                        int x1 = line.getQ().getX();
                        int x2 = line.getR().getX();
                        int y1 = line.getQ().getY();
                        int y2 = line.getR().getY();


                        codeOut = line.getR().bitcode;

                        // Find intersection point

                        //If line intersects top of window
                        if (Math.round(codeOut / top) == 1) {
                            try {
                                x = x1 + (x2 - x1) * (window.top - y1) / (y2 - y1);
                            } catch (ArithmeticException arithmeticException){
                                x = x1;
                            }
                            y = window.top;
                        }
                        //If line intersects bottom of window
                        else if (Math.round(codeOut / bottom) == 1) {
                            try{
                                x = x1 + (x2 - x1) * (window.bottom - y1) / (y2 - y1);
                            } catch (ArithmeticException arithmeticException){
                                x = x1;
                            }
                            y = window.bottom;

                        }
                        //If line intersects right of window
                        else if (Math.round(codeOut / right) == 1) {
                            x = window.right;
                            try {
                                y = y1 + (y2 - y1) * (window.right - x1) / (x2 - x1);
                            } catch (ArithmeticException arithmeticException){
                                y =y1;
                            }
                        }
                        //If line intersects left of window
                        else if (codeOut / left == 1) {
                            x = window.left;
                            try {
                                y = y1 + (y2 - y1) * (window.left - x1) / (x2 - x1);
                            } catch (ArithmeticException arithmeticException){
                                y = y1;
                            }
                        }
                        Point p = new Point(x, y);
                        assignBitCode(p);
                        line.setR(p);
                    }
                }
            }
            i++;
        }

        //remove out of bounds lines
        for(int j=0; j<outOfBounds.size(); j++){
            this.lines.remove(outOfBounds.get(j)-j);
        }
    }

    public void translation(int x, int y) {
        for (Line line : this.lines) {
            line.getQ().setX(line.getQ().getX() + x);
            line.getQ().setY(line.getQ().getY() + y);
            line.getR().setX(line.getR().getX() + x);
            line.getR().setY(line.getR().getY() + y);
        }
    }

    public void scale(double s) {
        for (Line line : this.lines) {
            line.getQ().setX((int) Math.round(line.getQ().getX() * s));
            line.getQ().setY((int) Math.round(line.getQ().getY() * s));
            line.getR().setX((int) Math.round(line.getR().getX() * s));
            line.getR().setY((int) Math.round(line.getR().getY() * s));
        }

    }

    public void rotation(int t) {
        double theta = Math.toRadians(t);

        for (Line line : this.lines) {
            double sint = Math.sin(theta);
            double cost = Math.cos(theta);

            int q_x = line.getQ().getX();
            int q_y = line.getQ().getY();
            int r_x = line.getR().getX();
            int r_y = line.getR().getY();

            //apply rotations
            int q_x_prime = (int) (q_x * cost - q_y*sint);
            int q_y_prime = (int) (q_x * sint + q_y*cost);
            int r_x_prime = (int) (r_x * cost - r_y*sint);
            int r_y_prime = (int) (r_x * sint + r_y*cost);

            //set rotations
            line.setQ(new Point(q_x_prime,q_y_prime));
            line.setR(new Point(r_x_prime,r_y_prime));

        }
    }

    public String writePbm() throws IOException {
        String pbmOutput = "P1\n";
        pbmOutput = pbmOutput.concat(String.format("# %s.pbm\n", this.fileName));
        pbmOutput = pbmOutput.concat(String.format("%d %d\n", this.pixels.length, this.pixels[0].length));
        for (int i = 0; i < length; i++) {
            int linecharcount = 0;
            for (int j = 0; j < width; j++) {
                if (this.pixels[j][i].isBlack()) {
                    pbmOutput = pbmOutput.concat("1 ");
                } else {
                    pbmOutput = pbmOutput.concat("0 ");
                }
                linecharcount++;
                linecharcount++;
                if (linecharcount % 70 == 0) {
                    pbmOutput = pbmOutput.concat("\n");
                }
            }
            pbmOutput = pbmOutput.concat("\n");
        }
        return pbmOutput;
    }

    public String writePs() {
        String psOutput = "/Line {moveto lineto stroke} bind def\n\n%%%BEGIN";
        //psOutput = psOutput.concat(String.valueOf(this.width));
        //psOutput = psOutput.concat(" ");
        //psOutput = psOutput = psOutput.concat(String.valueOf(this.length));
        //psOutput = psOutput.concat("] >> setpagedevice\n%%EndSetup\n\n%%%BEGIN\n");

        for (Line line : this.lines) {
            String output = String.valueOf(line.getQ().getX()) + " " + String.valueOf(line.getQ().getY()) + " " + String.valueOf(line.getR().getX()) + " " + String.valueOf(line.getR().getY()) + " Line\n";
            psOutput = psOutput.concat(output);
        }
        psOutput = psOutput.concat("%%%END");
        return psOutput;
    }

}