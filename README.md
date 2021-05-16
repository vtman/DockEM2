# DockEM2

<h1>Similarity of vectors</h1>

Suppose there be two 1D vectors of size <i>n</i>: <img src="https://latex.codecogs.com/svg.latex?\inline&space;f_i,&space;g_i,&space;i=0,\ldots,{n-1}" title="f_i, g_i, i=0,\ldots,{n-1}" />. The correlation coefficient is defined as
<img src="https://latex.codecogs.com/svg.latex?R_{fg}\equiv&space;\frac{\sum_{i=0}^{n-1}&space;(f_i&space;-&space;\bar&space;f)(g_i&space;-\bar&space;g)}{\sqrt{\sum_{i=0}^{n-1}(f_i-\bar&space;f)^2}&space;\sqrt{\sum_{i=0}^{n-1}(g_i-\bar&space;g)^2}}" title="R_{fg}\equiv \frac{\sum_{i=0}^{n-1} (f_i - \bar f)(g_i -\bar g)}{\sqrt{\sum_{i=0}^{n-1}(f_i-\bar f)^2} \sqrt{\sum_{i=0}^{n-1}(g_i-\bar g)^2}}" />

Symbols <img src="https://latex.codecogs.com/svg.latex?\inline&space;\bar&space;f" title="\bar f" /> and <img src="https://latex.codecogs.com/svg.latex?\inline&space;\bar&space;g" title="\bar g" /> are the mean values for the vectors, i.e. 
<img src="https://latex.codecogs.com/svg.latex?\inline&space;\bar&space;f\equiv&space;\frac{1}{n}\sum_{i=0}^{n-1}&space;f_i" title="\bar f\equiv \frac{1}{n}\sum_{i=0}^{n-1} f_i" />.

In case of an arbitrary vector <img src="https://latex.codecogs.com/svg.latex?\inline&space;\tilde&space;g_i" title="\tilde g_i" />, we may define vector such that

<img src="https://latex.codecogs.com/svg.latex?g_i&space;\equiv&space;\frac{\tilde&space;g_i&space;-&space;\frac{1}{n}\sum_{k=0}^{n-1}\tilde&space;g_k}{\sqrt{\sum_{p=0}^{n-1}\left(\tilde&space;g_p&space;-&space;\frac{1}{n}\sum_{k=0}^{n-1}\tilde&space;g_k\right)^2}}." title="g_i \equiv \frac{\tilde g_i - \frac{1}{n}\sum_{k=0}^{n-1}\tilde g_k}{\sqrt{\sum_{p=0}^{n-1}\left(\tilde g_p - \frac{1}{n}\sum_{k=0}^{n-1}\tilde g_k\right)^2}}." />

In this case <img src="https://latex.codecogs.com/svg.latex?\inline&space;\bar&space;g&space;=&space;0" title="\bar g = 0" /> and 
<img src="https://latex.codecogs.com/svg.latex?\inline&space;\sum_{i=0}^{n-1}&space;g_i^2&space;=&space;1" title="\sum_{i=0}^{n-1} g_i^2 = 1" />.
Then the normalised correlation coefficient can be written as

<img src="https://latex.codecogs.com/svg.latex?R_{fg}\equiv&space;\frac{\sum_{i=0}^{n-1}&space;f_i&space;g_i}{\sqrt{\sum_{i=0}^{n-1}&space;f_i^2-&space;n\bar&space;f^2}}" title="R_{fg}\equiv \frac{\sum_{i=0}^{n-1} f_i g_i}{\sqrt{\sum_{i=0}^{n-1} f_i^2- n\bar f^2}}"/>

The higher value of the correlation coefficient is, the better similarity between the given signals is.

We may scale the signal <i>f</i> linearly, i.e. use <i>F<sub>i</sub>=a f<sub>i</sub> +b</i>, where <i>a</i> and <i>b</i> are constants. The correlation coefficient is the same. Unfortunately, there are always some experimental and processing errors. As the result, in case of low signal-to-noise values the correlation coefficient may attain high values. In order to avoid this situation, we re-define the correlation coefficient as

<img src="https://latex.codecogs.com/svg.latex?R_{fg}\equiv&space;\left\{\begin{array}{ll}&space;\frac{\sum_{i=0}^{n-1}&space;f_i&space;g_i}{\sqrt{\sum_{i=0}^{n-1}&space;f_i^2-&space;n\bar&space;f^2}},&space;&&space;\mathrm{if\phantom{0}}&space;\sum_{i=0}^{n-1}&space;f_i^2-&space;n\bar&space;f^2&space;>&space;\sigma_0^2,\\&space;0,&space;&&space;\mathrm{otherwise.}&space;\end{array}&space;\right." title="R_{fg}\equiv \left\{\begin{array}{ll} \frac{\sum_{i=0}^{n-1} f_i g_i}{\sqrt{\sum_{i=0}^{n-1} f_i^2- n\bar f^2}}, & \mathrm{if\phantom{0}} \sum_{i=0}^{n-1} f_i^2- n\bar f^2 > \sigma_0^2,\\ 0, & \mathrm{otherwise.} \end{array} \right." />

<h1>Correlation maps</h1>


<h1>Command line</h1>



-target C:\Temp2\xDock\0420\in\emd_1997.map -search C:\Temp2\xDock\0420\in\EL-int-start6.pdb -outputFolder C:\Temp2\xDock\0420\out -prefix test -threads 6 -angularResolution 10.0 -angleForEstimation 15.0 -center 50,50,50 -boxSize 70,70,70  -signalLevel 0.5 -maskLevel 0.2 -maskRadius 4.0 -maskPixel 0.3 -resolution 6.0 -outputReduction 1 -systemFolder C:\Temp2\DockEM2\systemFolder -nRelativeInfo 25 -minCCpeak 0.45 -sphericalMask yes
