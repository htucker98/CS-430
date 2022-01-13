import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.io.FileOutputStream;
import java.util.Scanner;

public class Image {
    public File smf;
    public File pbm;
    public String fileName;
    public int j; //The next argument is an integer lower bound in the x dimension of the view port in screen coordinates (0)
    public int k; //integer in lower bound in the y dimension of the view port in screen coordinates (0)
    public int o; //integer in upper bound in the x dimension of the view port in screen coordinates (500)
    public int p; //integer in upper bound in the y dimension of the view port in screen coordinates (500)
    public float x; //floating point x of Projection Reference floatPoint (PRP) in VRC coordinates (0.0)
    public float y; //floating point y of Projection Reference floatPoint (PRP) in VRC coordinates (0.0)
    public float z; //floating point z of Projection Reference floatPoint (PRP) in VRC coordinates (1.0)
    public float X; //floating point x of View Reference floatPoint (VRP) in world coordinates (0.0)
    public float Y; //floating point y of View Reference floatPoint (VRP) in world coordinates (0.0)
    public float Z; //floating point z of View Reference floatPoint (VRP) in world coordinates (0.0)//
    public float q; //floating point x of View Plane Normal vector (VPN) in world coordinates (0.0)
    public float r; //floating point y of View Plane Normal vector (VPN) in world coordinates (0.0)
    public float w;
    public float Q;
    public float R;
    public float W;
    public float u;
    public float v;
    public float U;
    public float V;
    public float F = 0.6f;
    public float B = -0.6f;
    public boolean P = false;
    public FloatWindow viewport;
    public FloatWindow world;
    public Window window;
    public int width = 501;
    public int length = 501;
    public Pixel[][] pixels = new Pixel[length][width];
    public ArrayList<Vertex> vertices = new ArrayList<>();
    public ArrayList<float[][]> homogenousCoordinates = new ArrayList<>();
    public ArrayList<float[][]> transformedCoordinates = new ArrayList<>();
    public ArrayList<Face> faces = new ArrayList<>();
    public float[][] normalizingMatrix;
    public ArrayList<Vertex> projectionVertices = new ArrayList<>();
    public ArrayList<floatPoint> floatPoints = new ArrayList<>();
    public ArrayList<Point> points = new ArrayList<Point>();
    public ArrayList<ArrayList<Line>> polygonLines = new ArrayList<>();

    public Image(String[] args) throws FileNotFoundException{
        j = 0; k=0; o=500; p=500;
        x= 0.0f; y= 0.0f; z= 1.0f;
        X= 0.0f; Y= 0.0f; Z= 0.0f;
        q= 0.0f; r= 0.0f; w= -1.0f;
        Q= 0.0f; R= 1.0f; W= 0.0f;
        u= -0.7f; v= -0.7f;
        U= 0.7f; V= 0.7f;
        String f = System.getProperty("user.dir") + "/" + "bound-lo-sphere.smf";
        String fileName = "bound-lo-sphere";

        // check if there is specific commandline arguments
        for(int i=0; i< args.length; i++){
            if(args[i].equals("-f")){
                f = System.getProperty("user.dir") + "/" + args[i+1] ;
                fileName = f.replace(".smf","");
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
            if(args[i].equals("-x")){
                x = Float.parseFloat(args[i+1]);
            }
            if(args[i].equals("-y")){
                y = Float.parseFloat(args[i+1]);
            }
            if(args[i].equals("-z")){
                z = Float.parseFloat(args[i+1]);
            }
            if(args[i].equals("-X")){
                X = Float.parseFloat(args[i+1]);
            }
            if(args[i].equals("-Y")){
                Y = Float.parseFloat(args[i+1]);
            }
            if(args[i].equals("-Z")){
                Z = Float.parseFloat(args[i+1]);
            }
            if(args[i].equals("-q")){
                q = Float.parseFloat(args[i+1]);
            }
            if(args[i].equals("-r")){
                r = Float.parseFloat(args[i+1]);
            }
            if(args[i].equals("-w")){
                w = Float.parseFloat(args[i+1]);
            }
            if(args[i].equals("-Q")){
                Q = Float.parseFloat(args[i+1]);
            }
            if(args[i].equals("-R")){
                R = Float.parseFloat(args[i+1]);
            }
            if(args[i].equals("-W")){
                W = Float.parseFloat(args[i+1]);
            }
            if(args[i].equals("-u")){
                u = Float.parseFloat(args[i+1]);
            }
            if(args[i].equals("-v")){
                v = Float.parseFloat(args[i+1]);
            }
            if(args[i].equals("-U")){
                U = Float.parseFloat(args[i+1]);
            }
            if(args[i].equals("-V")){
                V = Float.parseFloat(args[i+1]);
            }
            if(args[i].equals("-P")){
                P = true;
            }
            if(args[i].equals("-F")){
                F = Float.parseFloat(args[i+1]);
            }
            if(args[i].equals("-B")){
                B = Float.parseFloat(args[i+1]);
            }
        }

        viewport = new FloatWindow(j,k,o,p);
        window = new Window(j,k,o,p);
        //set world window
        if(P){
            world = new FloatWindow(-1,-1,1,1);
        }
        else {
            float d = (z) / (B - z);
            d = java.lang.Math.abs(d);
            world = new FloatWindow(-(d), -(d), d, d);
        }


        //check if SFM file exists
        try {
            this.smf = new File(f);
            if (!this.smf.exists()) {
                throw new FileNotFoundException();
            }
        }
        catch (FileNotFoundException e) {
            e.printStackTrace();
            return;
        }

        //set fileName
        this.fileName = fileName;

        //create output pbm file
        this.pbm = new File(this.fileName + ".pbm");
        this.pixels = new Pixel[width][length];
        for (int i = 0; i <(this.width); i++) {
            for (int j = 0; j <(this.length); j++) {
                pixels[i][j] = new Pixel();
            }
        }
    }

    public void setPoints(){
        for(Vertex vertex : projectionVertices){
            floatPoints.add(new floatPoint(vertex.x,vertex.y));
        }
    }

    public void floatPointToPoint(){
        for (floatPoint floatPoint : floatPoints){
            int x = (int) floatPoint.x;
            int y = (int) floatPoint.y;
            points.add(new Point(x,y));

        }
    }

    public void facesToLines(){
        for (Face face: faces){
            ArrayList<Line> lines = new ArrayList<>();
            lines.add(new Line(points.get(face.x),points.get(face.y)));
            lines.add(new Line(points.get(face.y),points.get(face.z)));
            lines.add(new Line(points.get(face.z),points.get(face.x)));
            this.polygonLines.add(lines);
        }
    }

    public float[][] dotProduct(float[][] A, float[][] v){
        float [][] M = new float[A.length][v[0].length];

        for (int i=0; i<A.length;i++){
            for(int j=0; j<v[i].length; j++) {
                float dot = 0;
                for (int k = 0; k < A[i].length; k++) {
                    dot += (A[i][k] * v[k][j]);
                }
                M[i][j] = dot;
            }
        }
        return M;
    }
    public float[][] crossProduct(float[][]a, float[][] b){
        float x = (a[1][0]*b[2][0])-(a[2][0]*b[1][0]);
        float y = (a[2][0]*b[0][0])-(a[0][0]*b[2][0]);
        float z = (a[0][0]*b[1][0])-(a[1][0]*b[0][0]);
        float[][] v = {{x},{y},{z},{1}};
        return v;
    }
    public float magnitude(float[][] v){
        return (float) Math.sqrt((v[0][0]*v[0][0])+(v[1][0]*v[1][0])+(v[2][0]*v[2][0]));
    }
    public void homogenousCoordinates(){
        for (Vertex vertex : this.vertices){
            float[][] homoCord =  new float[][]{{vertex.x},{vertex.y},{vertex.z},{1}};
            this.homogenousCoordinates.add(homoCord);
        }
    }

    public void writePixel(int x, int y) {
        pixels[x][y].setBlack();
    }

    public void Scan() throws Exception {

        Scanner scanner = new Scanner(this.smf);

        //read file line by line
        String line = "";
        while (scanner.hasNextLine()) {
            line = scanner.nextLine();

            //Scan for vertices
            if (line.contains("v")) {
                String[] floats = line.split(" ");
                float x = Float.parseFloat(floats[1]);
                float y = Float.parseFloat(floats[2]);
                float z = Float.parseFloat(floats[3]);

                vertices.add(new Vertex(x,y,z));
            }

            //Scan for faces
            if (line.contains("f")) {
                String[] ints = line.split(" ");
                int x = Integer.parseInt(ints[1]);
                int y = Integer.parseInt(ints[2]);
                int z = Integer.parseInt(ints[3]);

                faces.add(new Face(x,y,z));
            }
        }
    }

    public void normalizingTransformation(){
        //Translate to VRP to origin
        float T[][] = {{1,0,0,-(X)},{0,1,0,-(Y)},{0,0,1,-(Z)},{0,0,0,1}};

        //Rotate

        float z_Mag = magnitude(new float[][]{{q},{r},{w}});
        float z1= q/z_Mag; float z2=r/ z_Mag; float z3=w/ z_Mag;
        float[][] Rz = {{z1},{z2},{z3}};

        float[][] vupxRz = crossProduct(new float[][]{{Q},{R},{W}},Rz);
        float x_Mag = magnitude(vupxRz);
        float x1=vupxRz[0][0]/x_Mag; float x2=vupxRz[1][0]; float x3=vupxRz[2][0];
        float[][] Rx = {{x1},{x2},{x3}};

        float[][] Ry = crossProduct(Rz,Rx);

        float[][] R = {{Rx[0][0],Rx[1][0],Rx[2][0],0},{Ry[0][0],Ry[1][0],Ry[2][0],0},{Rz[0][0],Rz[1][0],Rz[2][0],0},{0,0,0,1}};

        float[][] N = dotProduct(R,T);

        //parallel projection
        if(P){
            //Shear
            float S[][] = {{1,0,((0.5f*(U+u))-x)/z,0},{0,1,((0.5f*(V+v))-y)/z,0},{0,0,1,0},{0,0,0,1}};

            N = dotProduct(S,N);

            //Translate center to origin
            float C[][] = {{1,0,0,-(U+u)/2},{0,1,0,-(V+v)/2},{0,0,1,-(F)},{0,0,0,1}};

            N= dotProduct(C,N);

            //Scale
            float Sc[][] = {{2/(U-u),0,0,0},{0,2/(V-v),0,0},{0,0,1/(F-B)},{0,0,0,1}};

            N= dotProduct(Sc,N);
            normalizingMatrix = N;


        }
        else{
            //Translate
            float T2[][] = {{1,0,0,-(x)},{0,1,0,-(y)},{0,0,1,-(z)},{0,0,0,1}};

            N = dotProduct(T2,N);

            //Shear
            float S[][] = {{1,0,((0.5f*(U+u))-x)/z,0},{0,1,((0.5f*(V+v))-y)/z,0},{0,0,1,0},{0,0,0,1}};

            N= dotProduct(S,N);

            //Scale
            float Sc[][] = {{(2*z)/((U-u)*(z-B)),0,0,0},{0,(2*z)/((V-v)*(z-B)),0,0},{0,0,1/(z-B)},{0,0,0,1}};

            N= dotProduct(Sc,N);
            normalizingMatrix = N;
        }
        //apply transformation with normalizing matrix
        for (float[][] hCord: this.homogenousCoordinates){
            float[][] result = this.dotProduct(normalizingMatrix,hCord);
            this.transformedCoordinates.add(result);
        }
    }

    public void projectionTransformation(){
        //referenced equations: https://www.mauriciopoppe.com/notes/computer-graphics/viewing/projection-transform/
        float d = (z)/(B-z);
        float[][] M;

        //perspetive projection
        if(!P) {
            //return new float[][]{{(2*F)/(world.right-world.left),0,(world.right+world.left)/(world.right-world.left),0},{0,(2*F)/(world.top-world.bottom),(world.top+world.bottom)/(world.top-world.bottom),0},{0,0,-(B+F)/(B-F),-(B*F)/(B-F)},{0,0,-1,0}};
            M = new float[][]{{1, 0, 0, 0}, {0, 1, 0, 0}, {0, 0, 1, 0}, {0, 0, 1 / d, 0}};
        }

        //parallel projection
        else{
            M =  new float[][]{{1,0,0,0},{0,1,0,0},{0,0,0,0},{0,0,0,1}};
        }
        for (float[][] coord : this.transformedCoordinates){
            float[][] projCord = this.dotProduct(M,coord);
            this.projectionVertices.add(new Vertex(projCord[0][0]/projCord[3][0],projCord[1][0]/projCord[3][0],projCord[2][0]/(projCord[3][0])));
        }
    }

    public void viewportTransformation(){
        this.setPoints();
        translation(-world.left,-world.bottom);
        scale((viewport.right-viewport.left)/(world.right-world.left), ((viewport.top-viewport.bottom)/(world.top-world.bottom)));
        translation(viewport.left,viewport.bottom);
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



    public void sutherlandHodgman(Window window) {
        this.floatPointToPoint();
        this.facesToLines();
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
        for (ArrayList<Line> polygon : this.polygonLines) {
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
        this.polygonLines = topClippedPolygons;

    }

    public Point intersection(Point v_i, Point v_i1, String clipBoundary){
        int x = 0;
        int y = 0;
        int x1 = v_i.x;
        int y1 = v_i.y;
        int x2 = v_i1.x;
        int y2 = v_i1.y;

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

    public void translation(float x, float y) {

        for (floatPoint point : this.floatPoints) {
            point.x = (point.x + x);
            point.y = (point.y + y);
        }
    }

    public void scale(float s1, float s2) {
        for (floatPoint point : this.floatPoints) {
            point.x = point.x * s1;
            point.y = point.y * s2;
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
        psOutput = psOutput.concat("] >> setpagedevice\n%%EndSetup\n\n/Line {moveto lineto stroke} bind def\n\n%%%BEGIN\n");

        for (ArrayList<Line> polygon : this.polygonLines) {
            for (Line line : polygon) {
                String output = String.valueOf(line.getQ().x) + " " + String.valueOf(line.getQ().y) + " " + String.valueOf(line.getR().x) + " " + String.valueOf(line.getR().y) + " Line\n";
                psOutput = psOutput.concat(output);
            }
        }
        psOutput = psOutput.concat("%%%END");
        return psOutput;
    }

}

