import java.io.*;
import java.lang.reflect.Array;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;


public class Image {
    public File ps;
    public BufferedWriter pbm;
    public String fileName;
    public String outputFileName = "out.pbm";
    public ArrayList<ArrayList<Point>> polygonsPoint;
    public ArrayList<ArrayList<Line>> polygonsLine = new ArrayList<ArrayList<Line>>();
    public ArrayList<Line> lines;
    public Pixel[][] pixels;
    public Window world;
    public Window viewport;
    public int width;
    public int length;
    public int viewportWidth;
    public int viewportLength;
    public int a; public int b; public int c; public int d; public int j; public int k; public int o; public int p;


    public Image(String[] args) throws IOException {
        String f = System.getProperty("user.dir") + "/" + "hw3_split.ps";
        a = 0;
        b = 0;
        c = 250;
        d = 250;
        j = 0;
        k = 0;
        o = 200;
        p = 200;
        String fileName = "hw3_split";

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
                b = b/2;
            }
            if(args[i].equals("-c")){
                c = Integer.parseInt(args[i+1]);
            }
            if(args[i].equals("-d")){
                d = Integer.parseInt(args[i+1]);

            }
            if(args[i].equals("-j")){
                j = Integer.parseInt(args[i+1]);
            }
            if(args[i].equals("-k")){
                k = Integer.parseInt(args[i+1]);
            }
            if(args[i].equals("-o")){
                o = Integer.parseInt(args[i+1]);
            }
            if(args[i].equals("-p")){
                p = Integer.parseInt(args[i+1]);
            }
            if(args[i].equals(">")){
                outputFileName = args[i+1];
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

        //create buffered writer for output file
        this.pbm = new BufferedWriter(new FileWriter(this.outputFileName,false));

        //set fileName
        this.fileName = fileName;

        //intialize class attributes
        this.polygonsPoint = new ArrayList<ArrayList<Point>>();
        this.lines = new ArrayList<Line>();
        this.width = 501;
        this.length = 501;
        this.viewportWidth = (o-j)+1;
        this.viewportLength = (p-k)+1;

        this.world = new Window(a, b, c, d);
        this.viewport = new Window(j,k,o,p);
        this.pixels = new Pixel[this.width][this.length];
        for (int h = 0; h <this.width; h++) {
            for (int i = 0; i <(this.length); i++) {
                pixels[h][i] = new Pixel();
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
        ArrayList<Point> points = new ArrayList<>();
        while (!line.contains("%%%END")) {
            line = scanner.readLine();

            //Scan for polygon points
            if (line.contains("moveto")) {
                String[] ints = line.split(" ");
                int x; int y;
                x = Integer.parseInt(ints[0]);
                y = Integer.parseInt(ints[1]);

                points.add(new Point(x, y));
            }

            if (line.contains("lineto")){
                String[] ints = line.split("\\s+");
                int x; int y;
                x = Integer.parseInt(ints[0]);
                y = Integer.parseInt(ints[1]);

                points.add(new Point(x, y));
            }

            // if end of polygon remove extra point
            if (line.contains("stroke")){
                points.remove(points.size()-1);
                //add list of points to polygon list
                polygonsPoint.add(points);
                //reset point list
                points = new ArrayList<>();
            }
        }
    }

    public void pointsToLine(){
        this.polygonsLine = new ArrayList<>();
        for (ArrayList<Point> polygon : this.polygonsPoint) {
            ArrayList<Line> lines = new ArrayList<>();
            for (int j = 0; j < polygon.size()-1; j++) {
                lines.add(new Line(polygon.get(j), polygon.get(j + 1)));
            }
            //final edge
            lines.add(new Line(polygon.get(polygon.size()-1), polygon.get(0)));
            this.polygonsLine.add(lines);
        }
        linesToPoints();
    }

    public void linesToPoints(){
        for (ArrayList<Line> polygon : this.polygonsLine) {
            ArrayList<Point> points = new ArrayList<>();
            points.add(polygon.get(0).getQ());
            for (int j = 0; j < polygon.size()-1; j++) {
                if(!points.get(points.size()-1).equals(polygon.get(j).getQ())){
                    points.add(polygon.get(j).getQ());
                }
                if(!points.get(points.size()-1).equals(polygon.get(j).getR())){
                    points.add(polygon.get(j).getR());
                }

            }
            this.polygonsPoint.add(points);
        }
    }


    public void bresenham() {
            // refrenced at; http://tech-algorithm.com/articles/drawing-line-using-bresenham-algorithm/
        for (ArrayList<Line> polygon : polygonsLine)
            for (Line line : polygon) {
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
                        try {
                            this.writePixel(Math.abs(x)-1, Math.abs(this.length - 1) - (y - (this.world.bottom)));
                        }
                        catch (Exception e){
                            // don't write out of bound pixel
                        }
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


    public void findNextIntersection(int index, ArrayList<Line> polygon, String clipBoundary, Window window){

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
                    //line completely out of bounds don't add

                }

                //Case 4 : Output intersection point & vi+1 to P’’
                else if (clipBoundary >= v_i.getX() && clipBoundary <= v_i1.getX()) {
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
                else if (clipBoundary >= v_i.getY() && clipBoundary <= v_i1.getY()) {
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
                else if (clipBoundary <= v_i.getX() && clipBoundary >= v_i1.getX()) {
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
                else if (clipBoundary <= v_i.getY() && clipBoundary >= v_i1.getY()) {
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



    public void worldToViewport(){
        //refrence: https://www.javatpoint.com/computer-graphics-window-to-viewport-co-ordinate-transformation
        //translate to world origin
        this.translation(-world.left, -world.bottom);

        //scale to viewport
        double x_scale = ((double)(viewport.right-viewport.left)) / ((double) (world.right-world.left));
        double y_scale = ((double)(viewport.top-viewport.bottom)) / ((double)(world.top-world.bottom));

        for (ArrayList<Point> polygon: this.polygonsPoint) {
            for (Point point : polygon) {
                point.setX((int) Math.round(point.getX() * x_scale));
                point.setY((int) Math.round(point.getY() * y_scale));
            }
        }

        //translate to viewport origin
        this.translation(viewport.left, viewport.bottom);


    }

    public void polygonEndCleanUp() {
        for (int i=0; i<this.polygonsLine.size()-1; i++) {
            if (this.polygonsLine.get(i).size() != 0) {
                for(int j=0; j<this.polygonsLine.get(i).size()-1; j++){

                    Point a = this.polygonsLine.get(i).get(j).getR();
                    Point b = this.polygonsLine.get(i).get(j+1).getQ();
                    if (!a.equals(b)) {
                        //check if x or y value is at least equals
                        if (a.getX() == b.getX() || a.getY() == b.getY()) {
                            this.polygonsLine.get(i).add(j + 1, new Line(a, b));
                        } else {
                            //top left corner
                            if(a.x>b.x && a.y>b.y){
                                this.polygonsLine.get(i).add(j+1, new Line(a,new Point(b.x,a.y)));
                                this.polygonsLine.get(i).add(j+2, new Line(new Point(b.x,a.y),b));
                            }

                            // bottom right corner
                            if(a.x<b.x && a.y<b.y){
                                this.polygonsLine.get(i).add(j+1, new Line(a,new Point(a.x,b.y)));
                                this.polygonsLine.get(i).add(j+2, new Line(new Point(a.x,b.y),b));
                            }
                        }
                    }

                }

                //check first and last lines
                Point a = this.polygonsLine.get(i).get(this.polygonsLine.get(i).size()-1).getR();
                Point b = this.polygonsLine.get(i).get(0).getQ();
                if (!a.equals(b)){
                    //check if x or y value is at least equals
                    if(a.getX() == b.getX() || a.getY() == b.getY()){
                        this.polygonsLine.get(i).add(this.polygonsLine.get(i).size(), new Line(a,b));
                    }
                    else {
                        //top left corner
                        if(a.x>b.x && a.y>b.y){
                            this.polygonsLine.get(i).add(this.polygonsLine.get(i).size(), new Line(a,new Point(b.x,a.y)));
                            this.polygonsLine.get(i).add(this.polygonsLine.get(i).size(), new Line(new Point(b.x,a.y),b));
                        }

                        // bottom right corner
                        if(a.x<b.x && a.y<b.y){
                            this.polygonsLine.get(i).add(this.polygonsLine.get(i).size(), new Line(a,new Point(a.x,b.y)));
                            this.polygonsLine.get(i).add(this.polygonsLine.get(i).size(), new Line(new Point(a.x,b.y),b));
                        }
                    }
                }


            }
        }
    }





    public void fill() {
        //apply algorithm to each polygon
        for(int i=0; i<this.polygonsLine.size(); i++) {
            //for each scan line
            for (int scanLine = this.viewport.bottom; scanLine < this.viewport.top; scanLine++) {
                //create list of edges that crosses scan-line
                ArrayList<Line> scanLineIntersects = new ArrayList<>();
                for (Line line : this.polygonsLine.get(i)) {
                    //set up ymin & ymax values
                    int ymin;
                    int ymax;
                    if (line.getQ().getY() > line.getR().getY()) {
                        ymax = line.getQ().getY();
                        ymin = line.getR().getY();
                    } else {
                        ymax = line.getR().getY();
                        ymin = line.getQ().getY();
                    }

                    //if line is horizontal ignore
                    if (line.getR().getY() == line.getQ().getY()) {
                    }
                    //if ymax on scan-line ignore
                    else if (ymax == scanLine) {
                    }
                    //If ymin <= y < ymax add edge to scan-line y‘s edge list
                    else if (ymin <= scanLine && scanLine < ymax) {
                        scanLineIntersects.add(line);
                    }
                }
                //go through list of intersection lines
                ArrayList<Integer> Intersections = new ArrayList<>();

                //Calculate intersections with edges on list
                for (Line line : scanLineIntersects) {
                    int x = intersection(line,scanLine);
                    if (this.j <= x ) {
                        Intersections.add(intersection(line, scanLine));
                    }
                }

                //Sort intersections in x
                Collections.sort(Intersections);

                if (Intersections.size() > 1) {
                    for (int j = 0; j < Intersections.size(); j++) {
                        int x = Intersections.get(j);
                        if(j ==0) {
                            while (x != Intersections.get(j + 1)) {
                                try {
                                    this.writePixel(Math.abs(x-1), Math.abs(this.length - 1) - (scanLine - (this.world.bottom)));
                                    x++;
                                }
                                catch (Exception e){
                                    x++;
                                }
                            }

                        }

                        else if (j != Intersections.size()-1) {
                            while (x != Intersections.get(j + 1)) {
                                try{
                                    this.writePixel(Math.abs(x-1), Math.abs(this.length - 1) - (scanLine - (this.world.bottom)));
                                    x++;
                                }
                                catch(Exception e){
                                    x++;
                                }
                            }
                            try {
                                this.writePixel(Math.abs(x-1), Math.abs(this.length - 1) - (scanLine - (this.world.bottom)));
                            }
                            catch(Exception e) {
                                //do nothing
                            }
                        }
                        j++;
                    }
                }
            }
        }
    }


    public Point boundaryIntersection(Point v_i, Point v_i1, String clipBoundary, Window window){
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


    public int intersection(Line line, int scanLine){
        return line.getQ().getX() + (line.getR().getX() - line.getQ().getX()) * (scanLine - line.getQ().getY()) / (line.getR().getY() - line.getQ().getY());
    }

    public void rotation(int t) {
        double theta = Math.toRadians(t);
        for (ArrayList<Point> polygon: this.polygonsPoint) {
            for (Point point : polygon) {
                double sint = Math.sin(theta);
                double cost = Math.cos(theta);

                int x = point.getX();
                int y = point.getY();

                //apply rotations
                int x_prime = (int) (x * cost - y * sint);
                int y_prime = (int) (x * sint + y * cost);

                //set rotations
                point.setX(x_prime);
                point.setY(y_prime);
            }
        }
    }

    public void translation(int x, int y) {

        for (ArrayList<Point> polygon: this.polygonsPoint) {
            for (Point point : polygon) {
                point.setX(point.getX() + x);
                point.setY(point.getY() + y);
            }
        }
    }

    public void scale(double s) {
        for (ArrayList<Point> polygon: this.polygonsPoint) {
            for (Point point : polygon) {
                point.setX((int) Math.round(point.getX() * s));
                point.setY((int) Math.round(point.getY() * s));
            }
        }
    }

    public void writePbm(){
        System.out.print("P1\n");
        System.out.printf("# %s.pbm\n", this.fileName);
        System.out.printf("%d %d\n", this.pixels.length, this.pixels[0].length);
        for (int i = 0; i < length; i++) {
            int linecharcount = 0;
            for (int j = 0; j < width; j++) {
                if (this.pixels[j][i].isBlack()) {
                    System.out.print("1 ");
                } else {
                    System.out.print("0 ");
                }
                linecharcount++;
                linecharcount++;
                if (linecharcount % 70 == 0) {
                    System.out.print("\n");
                }
            }
            System.out.print("\n");
        }
    }

}