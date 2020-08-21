# How to fit a Bézier curve to a set of data?
Eric Rehm, Takuvik Joint Laboratory, Université Laval / CNRS, Québec

## Background
This question was originally asked in a StackOverflow posting
http://stackoverflow.com/questions/6299019/how-can-i-fit-a-b%C3%A9zier-curve-to-a-set-of-data/32723564#32723564

## Introduction
Say you have a set of data points to which you would like to fit cubic Bézier curves in a piecewise fashion.
There's a nice solution dating from 1995, complete with MATLAB code, that does this with G1 continuity. I have reimplemented the code, updating it to modern MATLAB (R2015b).

The solution constructs a piecewise G<sup>1</sup> cubic Bézier curve from cubic curve segments which have as their initial endpoints, or knot points, some of the data points.  The parameters for the curve are:  the knot points, the angles of the tangent vectors at the knot points, and the distances from each knot point to the adjacent control points.  The algoritm is developed and three solution curves are proposed:  Globally Optimized Only (GOO), Segmentally Optimized Only (SOO), and Segmented then Globally Optmized (SGO).  The algorithm can make a guess at the intiial endpoints (knot points) of each curve segment.

>  Lane, Edward J. *Fitting Data Using Piecewise G1 Cubic Bézier Curves*


>  Thesis, NAVAL POSTGRADUATE SCHOOL MONTEREY CA, 1995


> http://www.dtic.mil/dtic/tr/fulltext/u2/a298091.pdf

To use this you must, at minimum, specify the number of end (knot) points *n*.  This can be selected by trial and error, remembering that each cubie Bézier segment can represent at most two points of inflection.  Optionall, you can specify the knot points themselves, which increases the reliability of a fit. The thesis shows some pretty tough examples, some of which are shown below.   Note that Lane's approach guarantees G<sup>1</sup> continuity (directions of adjacent tangent vectors are identical) between the cubic Bézier segments, i.e., smooth joints. However, there can be discontinuities in curvature (changes in direction of second derivative).

## Brief description of the Lane 1995 algorithm
The algorithm proceeds as follows in the demonstration code (`BezierFitDemo.m`), where *Q* is an [Mx2] array of data points.

1. Select the number of knot points *n*.

        n = 3;  % Starting number of knot points

1. Produce an initial guess (IG) curve (`iguess0.m`). The IG curve is constructed geometrically from the knot points and the tangent(s).  The IG curve is not optimized with respect to all of the data. 
```matlab
Qt = Q';
[IG, k] = iguess0(Qt, n);
```
2. The IG curve is passed to a segment-wise optmization routine (`segop.m`).  This minimizes the distance between the each cubic Bézier segment and the data points, i.e., a piecewise least-squares solution or segmentally-optimized curve (SOC).
```matlab
SOC = segop(k, Qt, IG);
```
3. The SOC passed to a global optimizer (`globop.m`) which adjusts the control points for each cubic Bézier segment to minimize the total least-square error, producing the segmentally then globally optmized curve (GOC).
```matlab
GOC = globop(SOC, Qt, 0, k);
```
4.  The GOC is plotted using an internal plotting routine (`poplt.m`), along with the original data points.
```matlab
hgoc = poplt(GOC, Qt);
```
5.  The GOC is plotted again using the open source MATLAB [geom2D Toolbox](https://www.mathworks.com/matlabcentral/fileexchange/7844-geom2d) Bézier plotting routine `drawBezierCurve.m`. This demonstrates how to successively plot each cubic Bézier segment from the list of SGO control points extracted using cpoints(GOC).
```matlab
Cnew = cpoints(GOC)';
for i = 1:3:length(Cnew)-3
    disp(i:i+3)
        drawBezierCurve(Cnew(i:i+3,:));      % Fitted cubic bezier segment
end
```
6. Note that a Globally Optmized Only (GOO) curve could be produced by sending the initial guess (IG) to the global optmizer, e.g.,
```matlab
GOC = globop(IG, Qt, 0, k);
```
    This does not tend to produce a good result, so it is omitted in the demonstration code (`BezierFitDemo.m`)

## Examples
### Fit to data sampled from an existing Bézier curve (*n* = 3)
In this example, we generate a cubic Bézier with four control points. Then we try to fit it with *n*=3 knot points, leading to *n*-1=2 cubic Bézier
sections or 3*(*n*-1)+1=7 control points in all.  Since we are starting with a sampled Bézier, we have simply approximated the (exact) single subdivision of a single Bézier curve.  So, first we define control points to generate data.  These control points define a Bézier curve described in:
> *Solving the Nearest Point-on-Curve Problem* and
> *A Bezier Curve-Based Root-Finder*,
> both by Philip J. Schneider


> In  *Graphics Gems*, Academic Press, 1990


> http://www.realtimerendering.com/resources/GraphicsGems/gems/NearestPoint.c

```matlab
        C =[0     0;
            1     2;
            3     3;
            4     2];
```	    

Each piecewise cubic segment is drawin in a different color, along with the convex hull of the new control points found by the algorithm.
![Continuity](https://gitlab.com/erehm/PiecewiseG1BezierFit/raw/master/images/Example1.png "Credit: Eric Rehm, Université Laval")

### Piecewise cubic fit to a Lissajous figure (*n* = 3)
Here's an example of using just three knot points (chosen automatically by the code) to fit two cubic Bézier segments to a Lissajous figure.  If you look carefully at the top (where the Lissajous figure starts and stops), you will see a mismatch.  An exercise for the reader would be to specify the *n*=3 knot points directly, with the first and last knot points coinciding. 

Each piecewise cubic segment is drawin in a different color, along with the convex hull of the new control points found by the algorithm.
![Continuity](https://gitlab.com/erehm/PiecewiseG1BezierFit/raw/master/images/Example2.png "Credit: Eric Rehm, Université Laval")

### Piecewise cubic fit to two cycles of a sine wave (*n* = 5)
Another exercise for the reader would be to plot the RMS error in this approximation to a sine wave.  One could view this as a form of compression.

Each piecewise cubic segment is drawin in a different color, along with the convex hull of the new control points found by the algorithm. ![Continuity](https://gitlab.com/erehm/PiecewiseG1BezierFit/raw/master/images/Example3.png "Credit: Eric Rehm, Université Laval")

## Dependenciees
The demo code (`BezierFitDemo.m`) depends on the geom2D Toolbox to draw Bézier curves. geom2D requires MATLAB R2014b or later. 


https://www.mathworks.com/matlabcentral/fileexchange/7844-geom2d

## Background on continuity

* G<sup>0</sup>: Pieces are connected at endpoints.
* G<sup>1</sup>: Pieces are connected and have same unit tangent vector at endpoints.
* G<sup>2</sup>: Pieces are connected, have same unit tangent vector, and same curvature at endpoints.


G<sup>n</sup> implies all lower G<sup>i</sup>.


* C<sup>0</sup>: Pieces are connected at endpoints = G<sup>0</sup>.
* C<sup>1</sup>: Pieces are connected and have same unit velocity vector (tangent vector is not normalized in length).
* C<sup>2</sup>: Pieces are connected, have same unit velocity vector, and same acceleartion at endpoints.


C<sup>n</sup> implies all lower C<sup>i</sup>.


![Continuity](https://gitlab.com/erehm/PiecewiseG1BezierFit/raw/master/images/Continuity.jpg "Credit: Carlo Séquin, EECS, UC Berkeley")


