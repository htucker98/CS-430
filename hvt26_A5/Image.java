import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.io.FileOutputStream;
import java.util.Collections;
import java.util.Scanner;

public class Image {
    public File smfRed;
    public File smfGreen;
    public File smfBlue;
    public String pbm;
    public String fileNameRed;
    public String fileNameGreen;
    public String fileNameBlue;
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
    public ArrayList<Vertex> verticesRed = new ArrayList<>();
    public ArrayList<Vertex> verticesGreen = new ArrayList<>();
    public ArrayList<Vertex> verticesBlue = new ArrayList<>();
    public ArrayList<Face> facesRed = new ArrayList<>();
    public ArrayList<Face> facesGreen = new ArrayList<>();
    public ArrayList<Face> facesBlue = new ArrayList<>();
    public ArrayList<float[][]> homogenousCoordinatesRed = new ArrayList<>();
    public ArrayList<float[][]> homogenousCoordinatesGreen = new ArrayList<>();
    public ArrayList<float[][]> homogenousCoordinatesBlue = new ArrayList<>();
    public float[][] normalizingMatrix;
    public ArrayList<float[][]> transformedCoordinatesRed = new ArrayList<>();
    public ArrayList<float[][]> transformedCoordinatesGreen = new ArrayList<>();
    public ArrayList<float[][]> transformedCoordinatesBlue = new ArrayList<>();
    public ArrayList<Vertex> projectionVerticesRed = new ArrayList<>();
    public ArrayList<Vertex> projectionVerticesGreen = new ArrayList<>();
    public ArrayList<Vertex> projectionVerticesBlue = new ArrayList<>();
    public float[][] zBuffer;
    public ArrayList<floatPoint> floatPointsRed = new ArrayList<>();
    public ArrayList<floatPoint> floatPointsGreen = new ArrayList<>();
    public ArrayList<floatPoint> floatPointsBlue = new ArrayList<>();
    public ArrayList<Point> pointsRed = new ArrayList<Point>();
    public ArrayList<Point> pointsGreen = new ArrayList<Point>();
    public ArrayList<Point> pointsBlue = new ArrayList<Point>();
    public ArrayList<ArrayList<Line>> polygonLinesRed = new ArrayList<>();
    public ArrayList<ArrayList<Line>> polygonLinesGreen = new ArrayList<>();
    public ArrayList<ArrayList<Line>> polygonLinesBlue = new ArrayList<>();
    public int maxShade = 256;

    public Image(String[] args) throws FileNotFoundException{
        j = 0; k=0; o=500; p=500;
        x= 0.0f; y= 0.0f; z= 1.0f;
        X= 0.0f; Y= 0.0f; Z= 0.0f;
        q= 0.0f; r= 0.0f; w= -1.0f;
        Q= 0.0f; R= 1.0f; W= 0.0f;
        u= -0.7f; v= -0.7f;
        U= 0.7f; V= 0.7f;
        String fr = System.getProperty("user.dir") + "/" + "bound-sprellpsd.smf";
        String fileName = "bound-sprellpsd";

        // check if there is specific commandline arguments
        for(int i=0; i< args.length; i++){
            if(args[i].equals("-f")){
                fr = System.getProperty("user.dir") + "/" + args[i+1] ;
                fileNameRed = fr.replace(".smf","");
            }
            if(args[i].equals("-g")){
                String fg = System.getProperty("user.dir") + "/" + args[i+1] ;
                fileNameGreen = fg.replace(".smf","");
                smfGreen = new File(fg);
            }
            if(args[i].equals("-i")){
                String fb = System.getProperty("user.dir") + "/" + args[i+1] ;
                fileNameBlue = fb.replace(".smf","");
                smfBlue = new File(fb);
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
            if(args[i].equals(">")){
                pbm = (args[i+1]);
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


        smfRed = new File(fr);

        this.pixels = new Pixel[width][length];
        this.zBuffer = new float[width][length];
        for (int i = 0; i <(this.width); i++) {
            for (int j = 0; j <(this.length); j++) {
                pixels[i][j] = new Pixel(0,0,0);
            }
        }
    }

    public void setPoints(){
        for(Vertex vertex : projectionVerticesRed){
            floatPointsRed.add(new floatPoint(vertex.x,vertex.y));
        }
        if(projectionVerticesGreen != null){
            for(Vertex vertex : projectionVerticesGreen){
                floatPointsGreen.add(new floatPoint(vertex.x,vertex.y));
            }
        }
        if(projectionVerticesBlue != null){
            for(Vertex vertex : projectionVerticesBlue){
                floatPointsBlue.add(new floatPoint(vertex.x,vertex.y));
            }
        }
    }

    public void floatPointToPoint(){
        for (floatPoint floatPoint : floatPointsRed){
            int x = (int) floatPoint.x;
            int y = (int) floatPoint.y;
            pointsRed.add(new Point(x,y));
        }
        if(this.floatPointsGreen != null){
            for (floatPoint floatPoint : floatPointsGreen){
                int x = (int) floatPoint.x;
                int y = (int) floatPoint.y;
                pointsGreen.add(new Point(x,y));
            }
        }
        if(this.floatPointsBlue != null){
            for (floatPoint floatPoint : floatPointsBlue){
                int x = (int) floatPoint.x;
                int y = (int) floatPoint.y;
                pointsBlue.add(new Point(x,y));
            }
        }
    }

    public void facesToLines(){
        for (Face face: facesRed){
            ArrayList<Line> lines = new ArrayList<>();
            lines.add(new Line(pointsRed.get(face.x),pointsRed.get(face.y)));
            lines.add(new Line(pointsRed.get(face.y),pointsRed.get(face.z)));
            lines.add(new Line(pointsRed.get(face.z),pointsRed.get(face.x)));
            this.polygonLinesRed.add(lines);
        }
        if (facesGreen != null){
            for (Face face: facesGreen){
                ArrayList<Line> lines = new ArrayList<>();
                lines.add(new Line(pointsGreen.get(face.x),pointsGreen.get(face.y)));
                lines.add(new Line(pointsGreen.get(face.y),pointsGreen.get(face.z)));
                lines.add(new Line(pointsGreen.get(face.z),pointsGreen.get(face.x)));
                this.polygonLinesGreen.add(lines);
            }
        }
        if (facesBlue != null){
            for (Face face: facesBlue){
                ArrayList<Line> lines = new ArrayList<>();
                lines.add(new Line(pointsBlue.get(face.x),pointsBlue.get(face.y)));
                lines.add(new Line(pointsBlue.get(face.y),pointsBlue.get(face.z)));
                lines.add(new Line(pointsBlue.get(face.z),pointsBlue.get(face.x)));
                this.polygonLinesBlue.add(lines);
            }
        }
    }

    public int returnPixelShade(float zval){
        float factor = (this.F-zval)/(this.F-this.B);
        float val = 20f*factor;
        return (int) val;
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
        for (Vertex vertex : this.verticesRed){
            float[][] homoCord =  new float[][]{{vertex.x},{vertex.y},{vertex.z},{1}};
            this.homogenousCoordinatesRed.add(homoCord);
        }
        if(this.verticesGreen != null){
            for (Vertex vertex : this.verticesGreen){
                float[][] homoCord =  new float[][]{{vertex.x},{vertex.y},{vertex.z},{1}};
                this.homogenousCoordinatesGreen.add(homoCord);
            }
        }
        if(this.verticesBlue != null){
            for (Vertex vertex : this.verticesBlue){
                float[][] homoCord =  new float[][]{{vertex.x},{vertex.y},{vertex.z},{1}};
                this.homogenousCoordinatesBlue.add(homoCord);
            }
        }
    }

    public void Scan() throws Exception {

        Scanner scanner = new Scanner(this.smfRed);

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

                verticesRed.add(new Vertex(x,y,z));
            }

            //Scan for faces
            if (line.contains("f")) {
                String[] ints = line.split(" ");
                int x = Integer.parseInt(ints[1]);
                int y = Integer.parseInt(ints[2]);
                int z = Integer.parseInt(ints[3]);

                facesRed.add(new Face(x,y,z));
            }
        }

        if(this.smfGreen != null){
            scanner = new Scanner(this.smfGreen);
            //read file line by line
            line = "";
            while (scanner.hasNextLine()) {
                line = scanner.nextLine();

                //Scan for vertices
                if (line.contains("v")) {
                    String[] floats = line.split(" ");
                    float x = Float.parseFloat(floats[1]);
                    float y = Float.parseFloat(floats[2]);
                    float z = Float.parseFloat(floats[3]);

                    verticesGreen.add(new Vertex(x,y,z));
                }

                //Scan for faces
                if (line.contains("f")) {
                    String[] ints = line.split(" ");
                    int x = Integer.parseInt(ints[1]);
                    int y = Integer.parseInt(ints[2]);
                    int z = Integer.parseInt(ints[3]);

                    facesGreen.add(new Face(x,y,z));
                }
            }
        }

        if(this.smfBlue != null){
            scanner = new Scanner(this.smfBlue);
            //read file line by line
            line = "";
            while (scanner.hasNextLine()) {
                line = scanner.nextLine();

                //Scan for vertices
                if (line.contains("v")) {
                    String[] floats = line.split(" ");
                    float x = Float.parseFloat(floats[1]);
                    float y = Float.parseFloat(floats[2]);
                    float z = Float.parseFloat(floats[3]);

                    verticesBlue.add(new Vertex(x,y,z));
                }

                //Scan for faces
                if (line.contains("f")) {
                    String[] ints = line.split(" ");
                    int x = Integer.parseInt(ints[1]);
                    int y = Integer.parseInt(ints[2]);
                    int z = Integer.parseInt(ints[3]);

                    facesBlue.add(new Face(x,y,z));
                }
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
        for (float[][] hCord: this.homogenousCoordinatesRed){
            float[][] result = this.dotProduct(normalizingMatrix,hCord);
            this.transformedCoordinatesRed.add(result);
        }
        if(this.homogenousCoordinatesGreen != null){
            for (float[][] hCord: this.homogenousCoordinatesGreen){
                float[][] result = this.dotProduct(normalizingMatrix,hCord);
                this.transformedCoordinatesGreen.add(result);
            }
        }
        if(this.homogenousCoordinatesBlue != null){
            for (float[][] hCord: this.homogenousCoordinatesBlue){
                float[][] result = this.dotProduct(normalizingMatrix,hCord);
                this.transformedCoordinatesBlue.add(result);
            }
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
        for (float[][] coord : this.transformedCoordinatesRed){
            float[][] projCord = this.dotProduct(M,coord);
            this.projectionVerticesRed.add(new Vertex(projCord[0][0]/projCord[3][0],projCord[1][0]/projCord[3][0],projCord[2][0]/(projCord[3][0])));
        }
        if (this.transformedCoordinatesGreen != null) {
            for (float[][] coord : this.transformedCoordinatesGreen){
                float[][] projCord = this.dotProduct(M,coord);
                this.projectionVerticesGreen.add(new Vertex(projCord[0][0]/projCord[3][0],projCord[1][0]/projCord[3][0],projCord[2][0]/(projCord[3][0])));
            }
        }
        if (this.transformedCoordinatesBlue != null) {
            for (float[][] coord : this.transformedCoordinatesBlue){
                float[][] projCord = this.dotProduct(M,coord);
                this.projectionVerticesBlue.add(new Vertex(projCord[0][0]/projCord[3][0],projCord[1][0]/projCord[3][0],projCord[2][0]/(projCord[3][0])));
            }
        }
    }

    public void viewportTransformation(){
        translation(-world.left,-world.bottom);
        scale((viewport.right-viewport.left)/(world.right-world.left), ((viewport.top-viewport.bottom)/(world.top-world.bottom)));
        translation(viewport.left,viewport.bottom);
    }

    public void translation(float x, float y) {

        for (Vertex vertex : this.projectionVerticesRed) {
            vertex.x = (vertex.x + x);
            vertex.y = (vertex.y + y);
        }
        if(this.projectionVerticesGreen != null){
            for (Vertex vertex : this.projectionVerticesGreen) {
                vertex.x = (vertex.x + x);
                vertex.y = (vertex.y + y);
            }
        }
        if(this.projectionVerticesBlue != null){
            for (Vertex vertex : this.projectionVerticesBlue) {
                vertex.x = (vertex.x + x);
                vertex.y = (vertex.y + y);
            }
        }
    }

    public void scale(float s1, float s2) {
        for (Vertex vertex : this.projectionVerticesRed) {
            vertex.x = vertex.x * s1;
            vertex.y = vertex.y * s2;
        }
        if(this.projectionVerticesGreen != null){
            for (Vertex vertex : this.projectionVerticesGreen) {
                vertex.x = vertex.x * s1;
                vertex.y = vertex.y * s2;
            }
        }
        if(this.projectionVerticesBlue != null){
            for (Vertex vertex : this.projectionVerticesBlue) {
                vertex.x = vertex.x * s1;
                vertex.y = vertex.y * s2;
            }
        }
    }


    public void zbuf() {

        for (int y = 0; y < this.width; y++) {
            for (int x = 0; x < this.length; x++) {
                this.pixels[Math.abs((this.length-1)-(y-this.window.bottom))][x].setColor(0, 0, 0);
                this.zBuffer[Math.abs((this.length-1)-(y-this.window.bottom))][x] = this.B;
            }
        }
        //apply algorithm to each polygon
        for (Face face : this.facesRed) {

            Point p1, p2, p3;
            float z1, z2, z3, za, zb, zp;

            //if triangle is not on plan
            // e do not plot
            if (projectionVerticesRed.get(face.x).x == projectionVerticesRed.get(face.y).x && projectionVerticesRed.get(face.x).x == projectionVerticesRed.get(face.z).x || (projectionVerticesRed.get(face.x).y == projectionVerticesRed.get(face.y).y && projectionVerticesRed.get(face.x).y == projectionVerticesRed.get(face.z).y)) {

            } else if (projectionVerticesRed.get(face.x) == projectionVerticesRed.get(face.y) && projectionVerticesRed.get(face.x) == projectionVerticesRed.get(face.z)) {
                // TODO plot single vertex
            } else {

                p2 = null;
                z2 = 0f;
                p1 = null;
                z1 = 0f;
                p3 = null;
                z3 = 0f;
                // x left, y middle, z right x<=y y<=z
                if (projectionVerticesRed.get(face.x).x <= projectionVerticesRed.get(face.y).x && projectionVerticesRed.get(face.y).x <= projectionVerticesRed.get(face.z).x) {
                    p2 = new Point((int) projectionVerticesRed.get(face.x).x, (int) projectionVerticesRed.get(face.x).y);
                    z2 = transformedCoordinatesRed.get(face.x)[2][0];
                    p1 = new Point((int) projectionVerticesRed.get(face.y).x, (int) projectionVerticesRed.get(face.y).y);
                    z1 = transformedCoordinatesRed.get(face.y)[2][0];
                    p3 = new Point((int) projectionVerticesRed.get(face.z).x, (int) projectionVerticesRed.get(face.z).y);
                    z3 = transformedCoordinatesRed.get(face.z)[2][0];
                }
                // x left, z middle, y right x<=y z<=y x<z
                if (projectionVerticesRed.get(face.x).x <= projectionVerticesRed.get(face.y).x && projectionVerticesRed.get(face.z).x <= projectionVerticesRed.get(face.y).x && projectionVerticesRed.get(face.x).x <= projectionVerticesRed.get(face.z).x) {
                    p2 = new Point((int) projectionVerticesRed.get(face.x).x, (int) projectionVerticesRed.get(face.x).y);
                    z2 = transformedCoordinatesRed.get(face.x)[2][0];
                    p1 = new Point((int) projectionVerticesRed.get(face.z).x, (int) projectionVerticesRed.get(face.z).y);
                    z1 = transformedCoordinatesRed.get(face.z)[2][0];
                    p3 = new Point((int) projectionVerticesRed.get(face.y).x, (int) projectionVerticesRed.get(face.y).y);
                    z3 = transformedCoordinatesRed.get(face.y)[2][0];
                }
                // y left, x middle, z right y<=x x<=z
                if (projectionVerticesRed.get(face.y).x <= projectionVerticesRed.get(face.x).x && projectionVerticesRed.get(face.x).x <= projectionVerticesRed.get(face.z).x) {
                    p2 = new Point((int) projectionVerticesRed.get(face.y).x, (int) projectionVerticesRed.get(face.y).y);
                    z2 = transformedCoordinatesRed.get(face.y)[2][0];
                    p1 = new Point((int) projectionVerticesRed.get(face.x).x, (int) projectionVerticesRed.get(face.x).y);
                    z1 = transformedCoordinatesRed.get(face.x)[2][0];
                    p3 = new Point((int) projectionVerticesRed.get(face.z).x, (int) projectionVerticesRed.get(face.z).y);
                    z3 = transformedCoordinatesRed.get(face.z)[2][0];
                }
                // y left, z middle, x right y<=x z<=x y<z
                if (projectionVerticesRed.get(face.y).x <= projectionVerticesRed.get(face.x).x && projectionVerticesRed.get(face.z).x <= projectionVerticesRed.get(face.x).x && projectionVerticesRed.get(face.y).x <= projectionVerticesRed.get(face.z).x) {
                    p2 = new Point((int) projectionVerticesRed.get(face.y).x, (int) projectionVerticesRed.get(face.y).y);
                    z2 = transformedCoordinatesRed.get(face.y)[2][0];
                    p1 = new Point((int) projectionVerticesRed.get(face.z).x, (int) projectionVerticesRed.get(face.z).y);
                    z1 = transformedCoordinatesRed.get(face.z)[2][0];
                    p3 = new Point((int) projectionVerticesRed.get(face.x).x, (int) projectionVerticesRed.get(face.x).y);
                    z3 = transformedCoordinatesRed.get(face.x)[2][0];
                }
                // z left, x middle, y right z<=x x<=y
                if (projectionVerticesRed.get(face.z).x <= projectionVerticesRed.get(face.x).x && projectionVerticesRed.get(face.x).x <= projectionVerticesRed.get(face.y).x) {
                    p2 = new Point((int) projectionVerticesRed.get(face.z).x, (int) projectionVerticesRed.get(face.z).y);
                    z2 = transformedCoordinatesRed.get(face.z)[2][0];
                    p1 = new Point((int) projectionVerticesRed.get(face.x).x, (int) projectionVerticesRed.get(face.x).y);
                    z1 = transformedCoordinatesRed.get(face.x)[2][0];
                    p3 = new Point((int) projectionVerticesRed.get(face.y).x, (int) projectionVerticesRed.get(face.y).y);
                    z3 = transformedCoordinatesRed.get(face.y)[2][0];
                }
                // z left, y middle, x right z<=x y<=x z<y
                if (projectionVerticesRed.get(face.z).x <= projectionVerticesRed.get(face.x).x && projectionVerticesRed.get(face.y).x <= projectionVerticesRed.get(face.x).x && projectionVerticesRed.get(face.z).x <= projectionVerticesRed.get(face.y).x) {
                    p2 = new Point((int) projectionVerticesRed.get(face.z).x, (int) projectionVerticesRed.get(face.z).y);
                    z2 = transformedCoordinatesRed.get(face.z)[2][0];
                    p1 = new Point((int) projectionVerticesRed.get(face.y).x, (int) projectionVerticesRed.get(face.y).y);
                    z1 = transformedCoordinatesRed.get(face.y)[2][0];
                    p3 = new Point((int) projectionVerticesRed.get(face.x).x, (int) projectionVerticesRed.get(face.x).y);
                    z3 = transformedCoordinatesRed.get(face.x)[2][0];
                }

                ArrayList<Line> polygon = new ArrayList<>();
                polygon.add(new Line(p1, p2));
                polygon.add(new Line(p2, p3));
                polygon.add(new Line(p3, p1));

                try{
                    for (int scanline = p3.getY(); scanline <= p1.getY(); scanline++) {
                        za = z1 - ((z1 - z2) * ((float) p1.y - scanline) / ((float) p1.y - (float) p2.y));
                        zb = z1 - ((z1 - z3) * ((float) p1.y - scanline) / ((float) p1.y - (float) p3.y));
                        //x values at which line is intersected
                        ArrayList<Integer> Intersections = new ArrayList<>();
                        for (Line line : polygon) {
                            int ymin = Math.min(Math.min(p1.y, p2.y), p3.y);
                            int ymax = Math.max(Math.max(p1.y, p2.y), p3.y);

                            //if line is horizontal
                            if (line.getR().getY() == line.getQ().getY()) {
                                Intersections.add(line.getR().getX());
                                Intersections.add(line.getQ().getX());
                            } else if (ymax == scanline) {
                                Intersections.add(ymax);
                            }
                            //If ymin <= y < ymax add edge to scan-line y‘s edge list
                            else if (ymin <= scanline && scanline < ymax) {
                                Intersections.add(intersection(line, scanline));
                            }
                        }

                        //Sort intersections in x
                        Collections.sort(Intersections);

                        int start = Intersections.get(0);
                        int end = Intersections.get(1);
                        int x = start;
                        while (x != end + 1) {
                            //get polygon z-value at pixel location
                            zp = zb - ((zb - za) * ((float) end - (float) x) / ((float) end - (float) start));
                            try {
                                if (zp < this.F && zp > this.zBuffer[Math.abs((this.length-1)-(scanline-this.window.bottom))][x]) {
                                    this.zBuffer[Math.abs(Math.abs((this.length-1)-(scanline-this.window.bottom)))][x] = zp;
                                    int pix = this.returnPixelShade(zp);
                                    this.pixels[Math.abs((this.length-1)-(scanline-this.window.bottom))][x] = new Pixel(pix, 0, 0);
                                }
                            } catch (Exception e) {

                            }

                            x++;
                        }

                    }
                }
                catch (Exception e){

                }

            }
        }

        //apply algorithm to each green shape
        for (Face face : this.facesGreen) {

            Point p1, p2, p3;
            float z1, z2, z3, za, zb, zp;

            //if triangle is not on plan
            // e do not plot
            if (projectionVerticesGreen.get(face.x).x == projectionVerticesGreen.get(face.y).x && projectionVerticesGreen.get(face.x).x == projectionVerticesGreen.get(face.z).x || (projectionVerticesGreen.get(face.x).y == projectionVerticesGreen.get(face.y).y && projectionVerticesGreen.get(face.x).y == projectionVerticesGreen.get(face.z).y)) {

            } else if (projectionVerticesGreen.get(face.x) == projectionVerticesGreen.get(face.y) && projectionVerticesGreen.get(face.x) == projectionVerticesGreen.get(face.z)) {
                // TODO plot single vertex
            } else {

                p2 = null;
                z2 = 0f;
                p1 = null;
                z1 = 0f;
                p3 = null;
                z3 = 0f;
                // x left, y middle, z right x<=y y<=z
                if (projectionVerticesGreen.get(face.x).x <= projectionVerticesGreen.get(face.y).x && projectionVerticesGreen.get(face.y).x <= projectionVerticesGreen.get(face.z).x) {
                    p2 = new Point((int) projectionVerticesGreen.get(face.x).x, (int) projectionVerticesGreen.get(face.x).y);
                    z2 = transformedCoordinatesGreen.get(face.x)[2][0];
                    p1 = new Point((int) projectionVerticesGreen.get(face.y).x, (int) projectionVerticesGreen.get(face.y).y);
                    z1 = transformedCoordinatesGreen.get(face.y)[2][0];
                    p3 = new Point((int) projectionVerticesGreen.get(face.z).x, (int) projectionVerticesGreen.get(face.z).y);
                    z3 = transformedCoordinatesGreen.get(face.z)[2][0];
                }
                // x left, z middle, y right x<=y z<=y x<z
                if (projectionVerticesGreen.get(face.x).x <= projectionVerticesGreen.get(face.y).x && projectionVerticesGreen.get(face.z).x <= projectionVerticesGreen.get(face.y).x && projectionVerticesGreen.get(face.x).x <= projectionVerticesGreen.get(face.z).x) {
                    p2 = new Point((int) projectionVerticesGreen.get(face.x).x, (int) projectionVerticesGreen.get(face.x).y);
                    z2 = transformedCoordinatesGreen.get(face.x)[2][0];
                    p1 = new Point((int) projectionVerticesGreen.get(face.z).x, (int) projectionVerticesGreen.get(face.z).y);
                    z1 = transformedCoordinatesGreen.get(face.z)[2][0];
                    p3 = new Point((int) projectionVerticesGreen.get(face.y).x, (int) projectionVerticesGreen.get(face.y).y);
                    z3 = transformedCoordinatesGreen.get(face.y)[2][0];
                }
                // y left, x middle, z right y<=x x<=z
                if (projectionVerticesGreen.get(face.y).x <= projectionVerticesGreen.get(face.x).x && projectionVerticesGreen.get(face.x).x <= projectionVerticesGreen.get(face.z).x) {
                    p2 = new Point((int) projectionVerticesGreen.get(face.y).x, (int) projectionVerticesGreen.get(face.y).y);
                    z2 = transformedCoordinatesGreen.get(face.y)[2][0];
                    p1 = new Point((int) projectionVerticesGreen.get(face.x).x, (int) projectionVerticesGreen.get(face.x).y);
                    z1 = transformedCoordinatesGreen.get(face.x)[2][0];
                    p3 = new Point((int) projectionVerticesGreen.get(face.z).x, (int) projectionVerticesGreen.get(face.z).y);
                    z3 = transformedCoordinatesGreen.get(face.z)[2][0];
                }
                // y left, z middle, x right y<=x z<=x y<z
                if (projectionVerticesGreen.get(face.y).x <= projectionVerticesGreen.get(face.x).x && projectionVerticesGreen.get(face.z).x <= projectionVerticesGreen.get(face.x).x && projectionVerticesGreen.get(face.y).x <= projectionVerticesGreen.get(face.z).x) {
                    p2 = new Point((int) projectionVerticesGreen.get(face.y).x, (int) projectionVerticesGreen.get(face.y).y);
                    z2 = transformedCoordinatesGreen.get(face.y)[2][0];
                    p1 = new Point((int) projectionVerticesGreen.get(face.z).x, (int) projectionVerticesGreen.get(face.z).y);
                    z1 = transformedCoordinatesGreen.get(face.z)[2][0];
                    p3 = new Point((int) projectionVerticesGreen.get(face.x).x, (int) projectionVerticesGreen.get(face.x).y);
                    z3 = transformedCoordinatesGreen.get(face.x)[2][0];
                }
                // z left, x middle, y right z<=x x<=y
                if (projectionVerticesGreen.get(face.z).x <= projectionVerticesGreen.get(face.x).x && projectionVerticesGreen.get(face.x).x <= projectionVerticesGreen.get(face.y).x) {
                    p2 = new Point((int) projectionVerticesGreen.get(face.z).x, (int) projectionVerticesGreen.get(face.z).y);
                    z2 = transformedCoordinatesGreen.get(face.z)[2][0];
                    p1 = new Point((int) projectionVerticesGreen.get(face.x).x, (int) projectionVerticesGreen.get(face.x).y);
                    z1 = transformedCoordinatesGreen.get(face.x)[2][0];
                    p3 = new Point((int) projectionVerticesGreen.get(face.y).x, (int) projectionVerticesGreen.get(face.y).y);
                    z3 = transformedCoordinatesGreen.get(face.y)[2][0];
                }
                // z left, y middle, x right z<=x y<=x z<y
                if (projectionVerticesGreen.get(face.z).x <= projectionVerticesGreen.get(face.x).x && projectionVerticesGreen.get(face.y).x <= projectionVerticesGreen.get(face.x).x && projectionVerticesGreen.get(face.z).x <= projectionVerticesGreen.get(face.y).x) {
                    p2 = new Point((int) projectionVerticesGreen.get(face.z).x, (int) projectionVerticesGreen.get(face.z).y);
                    z2 = transformedCoordinatesGreen.get(face.z)[2][0];
                    p1 = new Point((int) projectionVerticesGreen.get(face.y).x, (int) projectionVerticesGreen.get(face.y).y);
                    z1 = transformedCoordinatesGreen.get(face.y)[2][0];
                    p3 = new Point((int) projectionVerticesGreen.get(face.x).x, (int) projectionVerticesGreen.get(face.x).y);
                    z3 = transformedCoordinatesGreen.get(face.x)[2][0];
                }

                ArrayList<Line> polygon = new ArrayList<>();
                polygon.add(new Line(p1, p2));
                polygon.add(new Line(p2, p3));
                polygon.add(new Line(p3, p1));

                try{
                    for (int scanline = p3.getY(); scanline <= p1.getY(); scanline++) {
                        za = z1 - ((z1 - z2) * ((float) p1.y - scanline) / ((float) p1.y - (float) p2.y));
                        zb = z1 - ((z1 - z3) * ((float) p1.y - scanline) / ((float) p1.y - (float) p3.y));
                        //x values at which line is intersected
                        ArrayList<Integer> Intersections = new ArrayList<>();
                        for (Line line : polygon) {
                            int ymin = Math.min(Math.min(p1.y, p2.y), p3.y);
                            int ymax = Math.max(Math.max(p1.y, p2.y), p3.y);

                            //if line is horizontal
                            if (line.getR().getY() == line.getQ().getY()) {
                                Intersections.add(line.getR().getX());
                                Intersections.add(line.getQ().getX());
                            } else if (ymax == scanline) {
                                Intersections.add(ymax);
                            }
                            //If ymin <= y < ymax add edge to scan-line y‘s edge list
                            else if (ymin <= scanline && scanline < ymax) {
                                Intersections.add(intersection(line, scanline));
                            }
                        }

                        //Sort intersections in x
                        Collections.sort(Intersections);

                        int start = Intersections.get(0);
                        int end = Intersections.get(1);
                        int x = start;
                        while (x != end + 1) {
                            //get polygon z-value at pixel location
                            zp = zb - ((zb - za) * ((float) end - (float) x) / ((float) end - (float) start));
                            try {
                                if (zp < this.F && zp > this.zBuffer[Math.abs((this.length-1)-(scanline-this.window.bottom))][x]) {
                                    this.zBuffer[Math.abs(Math.abs((this.length-1)-(scanline-this.window.bottom)))][x]= zp;
                                    int pix = this.returnPixelShade(zp);
                                    this.pixels[Math.abs((this.length-1)-(scanline-this.window.bottom))][x] = new Pixel(0, pix, 0);
                                }
                            } catch (Exception e) {

                            }

                            x++;
                        }

                    }
                }
                catch (Exception e){

                }

            }
        }

        for (Face face : this.facesBlue) {

            Point p1, p2, p3;
            float z1, z2, z3, za, zb, zp;

            //if triangle is not on plan
            // e do not plot
            if (projectionVerticesBlue.get(face.x).x == projectionVerticesBlue.get(face.y).x && projectionVerticesBlue.get(face.x).x == projectionVerticesBlue.get(face.z).x || (projectionVerticesBlue.get(face.x).y == projectionVerticesBlue.get(face.y).y && projectionVerticesBlue.get(face.x).y == projectionVerticesBlue.get(face.z).y)) {

            } else if (projectionVerticesBlue.get(face.x) == projectionVerticesBlue.get(face.y) && projectionVerticesBlue.get(face.x) == projectionVerticesBlue.get(face.z)) {
                // TODO plot single vertex
            } else {

                p2 = null;
                z2 = 0f;
                p1 = null;
                z1 = 0f;
                p3 = null;
                z3 = 0f;
                // x left, y middle, z right x<=y y<=z
                if (projectionVerticesBlue.get(face.x).x <= projectionVerticesBlue.get(face.y).x && projectionVerticesBlue.get(face.y).x <= projectionVerticesBlue.get(face.z).x) {
                    p2 = new Point((int) projectionVerticesBlue.get(face.x).x, (int) projectionVerticesBlue.get(face.x).y);
                    z2 = transformedCoordinatesBlue.get(face.x)[2][0];
                    p1 = new Point((int) projectionVerticesBlue.get(face.y).x, (int) projectionVerticesBlue.get(face.y).y);
                    z1 = transformedCoordinatesBlue.get(face.y)[2][0];
                    p3 = new Point((int) projectionVerticesBlue.get(face.z).x, (int) projectionVerticesBlue.get(face.z).y);
                    z3 = transformedCoordinatesBlue.get(face.z)[2][0];
                }
                // x left, z middle, y right x<=y z<=y x<z
                if (projectionVerticesBlue.get(face.x).x <= projectionVerticesBlue.get(face.y).x && projectionVerticesBlue.get(face.z).x <= projectionVerticesBlue.get(face.y).x && projectionVerticesBlue.get(face.x).x <= projectionVerticesBlue.get(face.z).x) {
                    p2 = new Point((int) projectionVerticesBlue.get(face.x).x, (int) projectionVerticesBlue.get(face.x).y);
                    z2 = transformedCoordinatesBlue.get(face.x)[2][0];
                    p1 = new Point((int) projectionVerticesBlue.get(face.z).x, (int) projectionVerticesBlue.get(face.z).y);
                    z1 = transformedCoordinatesBlue.get(face.z)[2][0];
                    p3 = new Point((int) projectionVerticesBlue.get(face.y).x, (int) projectionVerticesBlue.get(face.y).y);
                    z3 = transformedCoordinatesBlue.get(face.y)[2][0];
                }
                // y left, x middle, z right y<=x x<=z
                if (projectionVerticesBlue.get(face.y).x <= projectionVerticesBlue.get(face.x).x && projectionVerticesBlue.get(face.x).x <= projectionVerticesBlue.get(face.z).x) {
                    p2 = new Point((int) projectionVerticesBlue.get(face.y).x, (int) projectionVerticesBlue.get(face.y).y);
                    z2 = transformedCoordinatesBlue.get(face.y)[2][0];
                    p1 = new Point((int) projectionVerticesBlue.get(face.x).x, (int) projectionVerticesBlue.get(face.x).y);
                    z1 = transformedCoordinatesBlue.get(face.x)[2][0];
                    p3 = new Point((int) projectionVerticesBlue.get(face.z).x, (int) projectionVerticesBlue.get(face.z).y);
                    z3 = transformedCoordinatesBlue.get(face.z)[2][0];
                }
                // y left, z middle, x right y<=x z<=x y<z
                if (projectionVerticesBlue.get(face.y).x <= projectionVerticesBlue.get(face.x).x && projectionVerticesBlue.get(face.z).x <= projectionVerticesBlue.get(face.x).x && projectionVerticesBlue.get(face.y).x <= projectionVerticesBlue.get(face.z).x) {
                    p2 = new Point((int) projectionVerticesBlue.get(face.y).x, (int) projectionVerticesBlue.get(face.y).y);
                    z2 = transformedCoordinatesBlue.get(face.y)[2][0];
                    p1 = new Point((int) projectionVerticesBlue.get(face.z).x, (int) projectionVerticesBlue.get(face.z).y);
                    z1 = transformedCoordinatesBlue.get(face.z)[2][0];
                    p3 = new Point((int) projectionVerticesBlue.get(face.x).x, (int) projectionVerticesBlue.get(face.x).y);
                    z3 = transformedCoordinatesBlue.get(face.x)[2][0];
                }
                // z left, x middle, y right z<=x x<=y
                if (projectionVerticesBlue.get(face.z).x <= projectionVerticesBlue.get(face.x).x && projectionVerticesBlue.get(face.x).x <= projectionVerticesBlue.get(face.y).x) {
                    p2 = new Point((int) projectionVerticesBlue.get(face.z).x, (int) projectionVerticesBlue.get(face.z).y);
                    z2 = transformedCoordinatesBlue.get(face.z)[2][0];
                    p1 = new Point((int) projectionVerticesBlue.get(face.x).x, (int) projectionVerticesBlue.get(face.x).y);
                    z1 = transformedCoordinatesBlue.get(face.x)[2][0];
                    p3 = new Point((int) projectionVerticesBlue.get(face.y).x, (int) projectionVerticesBlue.get(face.y).y);
                    z3 = transformedCoordinatesBlue.get(face.y)[2][0];
                }
                // z left, y middle, x right z<=x y<=x z<y
                if (projectionVerticesBlue.get(face.z).x <= projectionVerticesBlue.get(face.x).x && projectionVerticesBlue.get(face.y).x <= projectionVerticesBlue.get(face.x).x && projectionVerticesBlue.get(face.z).x <= projectionVerticesBlue.get(face.y).x) {
                    p2 = new Point((int) projectionVerticesBlue.get(face.z).x, (int) projectionVerticesBlue.get(face.z).y);
                    z2 = transformedCoordinatesBlue.get(face.z)[2][0];
                    p1 = new Point((int) projectionVerticesBlue.get(face.y).x, (int) projectionVerticesBlue.get(face.y).y);
                    z1 = transformedCoordinatesBlue.get(face.y)[2][0];
                    p3 = new Point((int) projectionVerticesBlue.get(face.x).x, (int) projectionVerticesBlue.get(face.x).y);
                    z3 = transformedCoordinatesBlue.get(face.x)[2][0];
                }

                ArrayList<Line> polygon = new ArrayList<>();
                polygon.add(new Line(p1, p2));
                polygon.add(new Line(p2, p3));
                polygon.add(new Line(p3, p1));

                try{
                    for (int scanline = p3.getY(); scanline <= p1.getY(); scanline++) {
                        za = z1 - ((z1 - z2) * ((float) p1.y - scanline) / ((float) p1.y - (float) p2.y));
                        zb = z1 - ((z1 - z3) * ((float) p1.y - scanline) / ((float) p1.y - (float) p3.y));
                        //x values at which line is intersected
                        ArrayList<Integer> Intersections = new ArrayList<>();
                        for (Line line : polygon) {
                            int ymin = Math.min(Math.min(p1.y, p2.y), p3.y);
                            int ymax = Math.max(Math.max(p1.y, p2.y), p3.y);

                            //if line is horizontal
                            if (line.getR().getY() == line.getQ().getY()) {
                                Intersections.add(line.getR().getX());
                                Intersections.add(line.getQ().getX());
                            } else if (ymax == scanline) {
                                Intersections.add(ymax);
                            }
                            //If ymin <= y < ymax add edge to scan-line y‘s edge list
                            else if (ymin <= scanline && scanline < ymax) {
                                Intersections.add(intersection(line, scanline));
                            }
                        }

                        //Sort intersections in x
                        Collections.sort(Intersections);

                        int start = Intersections.get(0);
                        int end = Intersections.get(1);
                        int x = start;
                        while (x != end + 1) {
                            //get polygon z-value at pixel location
                            zp = zb - ((zb - za) * ((float) end - (float) x) / ((float) end - (float) start));
                            try {
                                if (zp < this.F && zp > this.zBuffer[Math.abs((this.length-1)-(scanline-this.window.bottom))][x]) {
                                    this.zBuffer[Math.abs((this.length-1)-(scanline-this.window.bottom))][x] = zp;
                                    int pix = this.returnPixelShade(zp);
                                    this.pixels[Math.abs((this.length-1)-(scanline-this.window.bottom))][x] = new Pixel(0, pix, 0);
                                }
                            } catch (Exception e) {

                            }

                            x++;
                        }

                    }
                }
                catch (Exception e){

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

    public void writePbm(){
        System.out.print("P3\n");
        System.out.printf("# %s.pbm\n", this.pbm);
        System.out.printf("%d %d\n", this.pixels.length, this.pixels[0].length);
        System.out.print("20\n");
        for (int i = 0; i < length; i++) {
            int linecharcount = 0;
            for (int j = 0; j < width; j++) {
                //print red
                System.out.printf("%d ",this.pixels[i][j].red);
                linecharcount+=3;
                if (linecharcount % 70 == 0) {
                    System.out.print("\n");
                }
                //print green
                System.out.printf("%d ",this.pixels[i][j].green);
                linecharcount+=3;
                if (linecharcount % 70 == 0) {
                    System.out.print("\n");
                }
                //print blue
                System.out.printf("%d ",this.pixels[i][j].blue);
                linecharcount+=3;
                if (linecharcount % 70 == 0) {
                    System.out.print("\n");
                }
                System.out.print("\t");
            }
            System.out.print("\n");
        }
    }

}

