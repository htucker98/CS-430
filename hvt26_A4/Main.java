public class Main {

    public static void main(String[] args) throws Exception {

        //check formatting of command line arguments
        Image image = new Image(args);

        image.Scan();
        image.homogenousCoordinates();
        image.normalizingTransformation();
        image.projectionTransformation();
        image.viewportTransformation();
        image.sutherlandHodgman(image.window);
        System.out.print(image.writePs());

    }
}
