#include "MatLabPlot.h"
#include "common.h"
using namespace std;
using namespace graph;

#ifndef M_PI
#define M_PI  3.14159265358979323846
#endif

double normal_pdf(double x, double mean, double sigma);




//http://anderberg.me/2016/08/01/c-variadic-templates/




int main() {
    
    vector<double> foo {25, 15, 5, -5, -15};
    vector<double> bar{1, 2, 3, 4, 5};
    
    vector<double> foo1 {1, 2, 3, 4, 5};
    vector<int> bar1{1, 2, 3, 4, 8};
 
    // Prepare sample data
    unsigned nPoint = 51;
    vector<double> x;
    x = linspace(0, 8, nPoint);
    vector<double> pdf11(nPoint);
    vector<double> pdf22(nPoint);
    double mu1 = 1;
    double mu2 = 2;
    double sigma1 = 1;
    double sigma2 = 2;

    for (unsigned i = 0; i < nPoint; ++i) {
        pdf11[i] = normal_pdf(x[i], mu1, sigma1);
        pdf22[i] = normal_pdf(x[i], mu2, sigma2);
    }


    xlabel("x");
    ylabel("pdf_{{/Symbol m},{/Symbol s}} (x)"); // enhanced texts

    legend({"({/Symbol m}, {/Symbol s}) = (1, 1)", "({/Symbol m}, {/Symbol s}) = (2, 2)"});
    
    grid(true); // turn on grids

    graph::plot(x, pdf11, "Color", "r", "LineWidth", 10, x, pdf22, "Color", "b", "LineWidth", 1);
    
    
    plot(foo, bar);
    
    plot(foo);
    
    return 0;
    


    vector<double> z = {3, 4, 5, 6, 9};
    vector<double> t = {1, 2, 3, 4, 5};

  

    return 0;
}



double normal_pdf(double x, double mean, double sigma) {
    x -= mean;
    x /= sigma;
    return 1.0 / sqrt(2 * M_PI) / sigma * exp(-x * x / 2);
}
