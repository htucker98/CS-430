
public class Main {

    public static void main(String[] args) throws Exception {

        //check formatting of command line arguments
        Image image = new Image(args);

        double s = 1.0;
        int m = 0;
        int n = 0;
        int r = 0;

        for(int i=0; i< args.length; i++){
            if(args[i].equals("-s")){
                s = Double.parseDouble(args[i+1]);
            }
            if(args[i].equals("-m")){
                m = Integer.parseInt(args[i+1]);
            }
            if(args[i].equals("-n")){
                n = Integer.parseInt(args[i+1]);
            }
            if(args[i].equals("-r")){
                r = Integer.parseInt(args[i+1]);
            }
        }

        image.Scan();
        if (s != 1.0){ image.scale(s);}
        if (r != 0){ image.rotation(r);}
        if (m != 0 || n!=0){ image.translation(m,n);}
        image.sutherlandHodgman();
        image.bresenham();
        System.out.print(image.writePs());



    }
}
