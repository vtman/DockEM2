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
Let there be two maps: <i>T</i> (target) and <i>S</i> (search). In cryo-EM the target map is often a density distribution obtained after reconstruction has been done. The search map may often be a simulated map, i.e. can be generated from atomic structures (PDB files). We aim to position the search map <i>S</i> onto the taret map <i>T</i>, so the corresponding common parts have the highest scores.  

For simplicity we consider a 1D case. The search signal may be defined on multiple segments. Thus we may extend both signals, so they are deived on the same segments. We may also introduce a mask vector <i>M</i> such that <i>M<sub>i</sub>=1</i> if <i>S<sub>i</sub></i> is known and set it zero otherwise. Note that <i>S<sub>i</sub> = S<sub>i</sub> M<sub>i</sub></i>. The total number of known values of the search vector is denoted as We define the number of known elements as <img src="https://latex.codecogs.com/svg.latex?\inline&space;w\equiv&space;\sum_{i=0}^{N-1}&space;M_i" title="w\equiv \sum_{i=0}^{N-1} M_i" />. The correlation coefficient can be written as

<img src="https://latex.codecogs.com/svg.latex?R^q_{TS}&space;=&space;\frac{\sum_{j=0}^{N-1}&space;T_{j&space;&plus;&space;q}&space;S_j&space;M_j}{&space;\sqrt{\sum_{j=0}^{N-1}&space;T^2_{j&plus;q}&space;M_j&space;-&space;w&space;\left(\frac{1}{w}\sum_{j=0}^{N-1}&space;T_{j&plus;q}M_j\right)^2}}." title="R^q_{TS} = \frac{\sum_{j=0}^{N-1} T_{j + q} S_j M_j}{ \sqrt{\sum_{j=0}^{N-1} T^2_{j+q} M_j - w \left(\frac{1}{w}\sum_{j=0}^{N-1} T_{j+q}M_j\right)^2}}." />

In case of two discrete real valued functions <img src="https://latex.codecogs.com/svg.latex?\inline&space;f_i,&space;g_i,&space;i=-\infty,\ldots,\infty" title="f_i, g_i, i=-\infty,\ldots,\infty" />  their cross-correlation is defined as

<img src="https://latex.codecogs.com/svg.latex?(f\star&space;g)[j]\equiv&space;\sum_{i=-\infty}^{\infty}&space;f_i&space;g_{i&space;&plus;&space;j}" title="(f\star g)[j]\equiv \sum_{i=-\infty}^{\infty} f_i g_{i + j}" />

So, we may rewrite the above formula as

<img src="https://latex.codecogs.com/svg.latex?R^q_{TS}&space;=&space;\frac{(S\star&space;T)[q]}{&space;\sqrt{(M\star&space;T^2)[q]&space;-&space;\frac{1}{w}&space;\left((M\star&space;T)[q]\right)^2}}." title="R^q_{TS} = \frac{(S\star T)[q]}{ \sqrt{(M\star T^2)[q] - \frac{1}{w} \left((M\star T)[q]\right)^2}}." />

<h1>Command line</h1>



-target C:\Temp2\xDock\0420\in\emd_1997.map -search C:\Temp2\xDock\0420\in\EL-int-start6.pdb -outputFolder C:\Temp2\xDock\0420\out -prefix test -threads 6 -angularResolution 10.0 -angleForEstimation 15.0 -center 50,50,50 -boxSize 70,70,70  -signalLevel 0.5 -maskLevel 0.2 -maskRadius 4.0 -maskPixel 0.3 -resolution 6.0 -outputReduction 1 -systemFolder C:\Temp2\DockEM2\systemFolder -nRelativeInfo 25 -minCCpeak 0.45 -sphericalMask yes
