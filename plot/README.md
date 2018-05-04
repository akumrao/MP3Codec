Work is on progress:

This is most simple and smaller Mpeg2 Layer 3 encoder. It uses hyprid MPEG Psychoacoustic Models, it uses both PAM1 and PAM2 features.

It compiles in all Linux, unix and embeeded OS having C++11



*******************************************************************************************************************************************
MPEG Psychoacoustic Models1. Matlab and Octave code
MPEG Psychoacoustic Models2. Matlab and Octave code

*******************************************************************************************************************************************


Graph plot is library to plot and analyze the MPEG Psychoacoustic Models
=======

gnuplot4.6 should be installed. PDF files will be generated 

A C++ library that allows you to plot just as in MATLAB 


Platforms
---------

Windows, Linux, OS X
C++11 
gnuplot 4.6 and above in `$PATH`


examples
--------------

### 1. Simple plot


graph::Plot curvePlot; 

curvePlot.plot({t,x});
curvePlot.exec();


### 2. Simple plot with legends


graph::Plot curvePlot; 

curvePlot.xlabel("x");
curvePlot.ylabel("pdf_{{/Symbol m},{/Symbol s}} (x)");  // enhanced texts

curvePlot.plot({ t,x1, t,x2 });  // plot two curves by cascading data pairs

curvePlot.legend({"{/Symbol l}=1",
                  "{/Symbol l}=4"});
curvePlot.grid(true);            // turn on grids

curvePlot.exec();
```


Texts can be formatted with subscripts (`"_"`), superscript (`"^"`), and breaklines (`"\\n"`).
Greek letters and other symbols are represented in different syntaxes from Latex. 


### 2. Simple plot with marker and color of graph

graph::Plot curvePlot; 
curvePlot.plot({ t,x1, t,x2 });



// multiple property setup in one statement
curvePlot.linespec(1, {{MarkerSize, "0.5"}, {Marker, "*"}});

// single property setup in multiple statements
curvePlot.linespec(2, Color, "b");
curvePlot.linespec(2, LineWidth, 2);
curvePlot.linespec(2, Marker, "none");

curvePlot.exec();

