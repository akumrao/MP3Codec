
A C++ library that allows you to plot just as in MATLAB 
=======

gnuplot4.6 should be installed. PDF files will be generated 




Platforms
---------

Windows, Linux, OS X
C++14 
gnuplot 4.6 and above in `$PATH`


examples
--------------

### 1. Simple plot


graph::Plot curvePlot; 

vector<double> x{1, 2, 3, 4, 5};
vector<double> y {25, 15, 5, -5, -15};
plot(x, y, LineStyle, "--" );



### 2. Simple plot with legends

label("x");
ylabel("y"); // enhanced texts

legend({"({/Symbol m}, {/Symbol s}) = (1, 1)", "({/Symbol m}, {/Symbol s}) = (2, 2)"});
grid(true); 
plot(x, y, "Color", "r", "LineWidth", 10, x, y2, "Color", "b", "LineWidth", 1);
   
