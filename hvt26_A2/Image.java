import java.io.File;
import java.io.FileReader;
import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;


public class Image {
    public File ps;
    public File pbm;
    public String fileName;
    public ArrayList<Point> points;
    public ArrayList<Line> lines;
    public Pixel[][] pixels;
    public Window window;
    public int width;
    public int length;


    public Image(String[] args) throws FileNotFoundException{
        String f = System.getProperty("user.dir") + "/" + "hw2_a.ps";
        int a = 0;
        int b = 0;
        int c = 499;
        int d = 499;
        String fileName = "hw2_a";

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
        this.points = new ArrayList<Point>();
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

            //Scan for polygon points
            if (line.contains("moveto")) {
                String[] ints = line.split(" ");
                int x = Integer.parseInt(ints[0]);
                int y = Integer.parseInt(ints[1]);
                points.add(new Point(x, y));
            }

            if (line.contains("lineto")){
                String[] ints = line.split(" ");
                int x = Integer.parseInt(ints[0]);
                int y = Integer.parseInt(ints[1]);
                points.add(new Point(x, y));
                points.add(new Point(x, y));
            }

            // if end of polygon remove extra point
            if (line.contains("stroke")){
                points.remove(points.size()-1);
            }
        }
    }

    public void pointsToLine(){
        for (int i = 0; i < this.points.size()-1; i = i + 2){
            this.lines.add(new Line(this.points.get(i),this.points.get(i+1)));
        }
    }

    public void bresenham() {
        pointsToLine();


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

    public void sutherlandHodgman(Window window) {
        //polygon clipping algorithm

        //define method variables
        boolean inside_vi;
        boolean inside_vi1;
        int clipBoundary;
        Point v_i;
        Point v_i1;


        ArrayList<ArrayList<Line>> leftClippedPolygons = new ArrayList<ArrayList<Line>>();
        //check points for LEFT clip boundary
        clipBoundary = window.left;
        for (ArrayList<Line> polygon : this.polygonsLine) {
            ArrayList<Line> leftClippedLines = new ArrayList<>();
            for (int j = 0; j < polygon.size(); j++) {

                v_i = polygon.get(j).getQ();
                v_i1 = polygon.get(j).getR();

                //Case 1 : Output vi+1 to P’’
                if (clipBoundary <= v_i.getX() && clipBoundary <= v_i1.getX()) {
                    //no change needed both points inside clip boundary
                    leftClippedLines.add(new Line(v_i,v_i1));
                }

                //Case 2 : Output intersection point to P’’
                else if (clipBoundary <= v_i.getX() && clipBoundary > v_i1.getX()) {
                    //find intersection with top boundary
                    v_i1 = boundaryIntersection(v_i, v_i1, "left", window);
                    leftClippedLines.add(new Line(v_i,v_i1));
                }

                //Case 3 : No output
                else if (clipBoundary > v_i.getX() && clipBoundary > v_i1.getX()) {
                    //line completely out of bounds do not add
                }

                //Case 4 : Output intersection point & vi+1 to P’’
                else if (clipBoundary > v_i.getX() && clipBoundary <= v_i1.getX()) {
                    //find intersection with top boundary
                    v_i = boundaryIntersection(v_i, v_i1, "left", window);
                    leftClippedLines.add(new Line(v_i,v_i1));
                }
            }
            leftClippedPolygons.add(leftClippedLines);
        }

        //check points for BOTTOM clip boundary
        ArrayList<ArrayList<Line>> bottomClippedPolygons = new ArrayList<ArrayList<Line>>();
        clipBoundary = window.bottom;
        for (ArrayList<Line> polygon : leftClippedPolygons) {
            ArrayList<Line> bottomClippedLines = new ArrayList<>();
            for (int j = 0; j < polygon.size(); j++) {

                v_i = polygon.get(j).getQ();
                v_i1 = polygon.get(j).getR();

                //Case 1 : Output vi+1 to P’’
                if (clipBoundary <= v_i.getY() && clipBoundary <= v_i1.getY()) {
                    //no change needed both points inside clip boundary
                    bottomClippedLines.add(new Line(v_i,v_i1));
                }

                //Case 2 : Output intersection point to P’’
                else if (clipBoundary <= v_i.getY() && clipBoundary > v_i1.getY()) {
                    //find intersection with top boundary
                    v_i1 = boundaryIntersection(v_i, v_i1, "bottom", window);
                    bottomClippedLines.add(new Line(v_i,v_i1));
                }

                //Case 3 : No output
                else if (clipBoundary > v_i.getY() && clipBoundary > v_i1.getY()) {
                }

                //Case 4 : Output intersection point & vi+1 to P’’
                else if (clipBoundary > v_i.getY() && clipBoundary <= v_i1.getY()) {
                    //find intersection with top boundary
                    v_i = boundaryIntersection(v_i, v_i1, "bottom", window);
                    bottomClippedLines.add(new Line(v_i,v_i1));
                }
            }
            bottomClippedPolygons.add(bottomClippedLines);
        }

        //check points for RIGHT clip boundary
        ArrayList<ArrayList<Line>> rightClippedPolygons = new ArrayList<ArrayList<Line>>();
        clipBoundary = window.right;
        for (ArrayList<Line> polygon : bottomClippedPolygons) {
            ArrayList<Line> rightClippedLines = new ArrayList<>();
            for (int j = 0; j < polygon.size(); j++) {

                v_i = polygon.get(j).getQ();
                v_i1 = polygon.get(j).getR();

                //Case 1 : Output vi+1 to P’’
                if (clipBoundary >= v_i.getX() && clipBoundary >= v_i1.getX()) {
                    //no change needed both points inside clip boundary
                    rightClippedLines.add(new Line(v_i,v_i1));
                }

                //Case 2 : Output intersection point to P’’
                else if (clipBoundary >= v_i.getX() && clipBoundary < v_i1.getX()) {
                    //find intersection with top boundary
                    v_i1 = boundaryIntersection(v_i, v_i1, "right", window);
                    rightClippedLines.add(new Line(v_i,v_i1));
                }

                //Case 3 : No output
                else if (clipBoundary < v_i.getX() && clipBoundary < v_i1.getX()) {
                }

                //Case 4 : Output intersection point & vi+1 to P’’
                else if (clipBoundary < v_i.getX() && clipBoundary >= v_i1.getX()) {
                    //find intersection with top boundary
                    v_i = boundaryIntersection(v_i, v_i1, "right", window);
                    rightClippedLines.add(new Line(v_i,v_i1));
                }
            }
            rightClippedPolygons.add(rightClippedLines);
        }

        //check points for TOP clip boundary
        ArrayList<ArrayList<Line>> topClippedPolygons = new ArrayList<ArrayList<Line>>();
        clipBoundary = window.top;
        for (ArrayList<Line> polygon : rightClippedPolygons) {
            ArrayList<Line> topClippedLines = new ArrayList<>();
            for (int j = 0; j < polygon.size(); j++) {

                v_i = polygon.get(j).getQ();
                v_i1 = polygon.get(j).getR();

                //Case 1 : Output vi+1 to P’’
                if (clipBoundary >= v_i.getY() && clipBoundary >= v_i1.getY()) {
                    //no change needed both points inside clip boundary
                    topClippedLines.add(new Line(v_i,v_i1));
                }

                //Case 2 : Output intersection point to P’’
                else if (clipBoundary >= v_i.getY() && clipBoundary < v_i1.getY()) {
                    //find intersection with top boundary
                    v_i1 = boundaryIntersection(v_i, v_i1, "top", window);
                    topClippedLines.add(new Line(v_i,v_i1));
                }

                //Case 3 : No output
                else if (clipBoundary < v_i.getY() && clipBoundary < v_i1.getY()) {

                }

                //Case 4 : Output intersection point & vi+1 to P’’
                else if (clipBoundary < v_i.getY() && clipBoundary >= v_i1.getY()) {
                    //find intersection with top boundary
                    v_i = boundaryIntersection(v_i, v_i1, "top", window);
                    topClippedLines.add(new Line(v_i,v_i1));
                }
            }
            topClippedPolygons.add(topClippedLines);
        }

        //set polygons to clipped polygons
        this.polygonsLine = topClippedPolygons;

    }

    public Point intersection(Point v_i, Point v_i1, String clipBoundary){
        int x = 0;
        int y = 0;
        int x1 = v_i.getX();
        int y1 = v_i.getY();
        int x2 = v_i1.getX();
        int y2 = v_i1.getY();

        // Find intersection point
        if (clipBoundary.equals("top")) {
            try {
                x = x1 + (x2 - x1) * (window.top - y1) / (y2 - y1);
            } catch (ArithmeticException ae) {
                x = x1;
            }

            y = window.top;
        }
        else if (clipBoundary.equals("bottom")) {
            try {
                x = x1 + (x2 - x1) * (window.bottom - y1) / (y2 - y1);
            } catch (ArithmeticException ae) {
                x = x1;
            }

            y = window.bottom;

        }
        else if (clipBoundary.equals("right")) {
            x = window.right;

            try {
                y = y1 + (y2 - y1) * (window.right - x1) / (x2 - x1);
            } catch (ArithmeticException ae) {
                y = y1;
            }

        }
        else if (clipBoundary.equals("left")){
            x = window.left;

            try {
                y = y1 + (y2 - y1) * (window.left - x1) / (x2 - x1);
            } catch (ArithmeticException ae) {
                y = y1;
            }

        }

        //return intersection point
        return new Point(x,y);
    }

    public void translation(int x, int y) {

        for (Point point : this.points) {
            point.setX(point.getX() + x);
            point.setY(point.getY() + y);
        }
    }

    public void scale(double s) {
        for (Point point : this.points) {
            point.setX((int) Math.round(point.getX() * s));
            point.setY((int) Math.round(point.getY() * s));
        }
    }

    public void rotation(int t) {
        double theta = Math.toRadians(t);


        for (Point point : this.points) {
            double sint = Math.sin(theta);
            double cost = Math.cos(theta);

            int x = point.getX();
            int y = point.getY();

            //apply rotations
            int x_prime = (int) (x * cost - y*sint);
            int y_prime = (int) (x * sint + y*cost);

            point.setX(x_prime); point.setY(y_prime);
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
        String psOutput = "%%BeginSetup\n<< /PageSize [";
        psOutput = psOutput.concat(String.valueOf(this.width));
        psOutput = psOutput.concat(" ");
        psOutput = psOutput.concat(String.valueOf(this.length));
        psOutput = psOutput.concat("] >> setpagedevice\n%%EndSetup\n/Line {moveto lineto stroke} bind def\n\n%%%BEGIN\n");

        for (Line line : this.lines) {
            String output = String.valueOf(line.getQ().getX()) + " " + String.valueOf(line.getQ().getY()) + " " + String.valueOf(line.getR().getX()) + " " + String.valueOf(line.getR().getY()) + " Line\n";
            psOutput = psOutput.concat(output);
        }
        psOutput = psOutput.concat("%%%END");
        return psOutput;
    }

}