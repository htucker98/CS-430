CS430 - A1 - Draw Clipped Lines

About: Program reads in a simplifed Postscript file (*.ps) that only contains lines and enables the user to perform transformations(translate, rotate, and scale). After, if any transformations are applied the image is clipped to the window dimensions and the bresenham algorithm is used to plot pixels according to a lines slope.

Language: Java 11
OS: macOS Mojave 10.14.16

Complile/Execution Notes:
The main() is located in Main.java.
To run simiply use makefile and run CG_hw1 as a shell script and indicate output with >
Example: ./CG_hw1 > out.ps

Line Clipping Location:
    File: Image.java
    Lines: 167-313