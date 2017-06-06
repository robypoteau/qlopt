import numpy as np

def format_data(v, noise):
    N = v.size
    odd = np.mod(N,2)
    v = np.sort(v)

    median = np.median(v)
    twentyfifth = np.percentile(v, 25)
    seventyfifth = np.percentile(v, 75)
    
    iqr = seventyfifth - twentyfifth
    upper_whisker = v[v<=seventyfifth+1.5*iqr].max()
    lower_whisker = v[v>=twentyfifth-1.5*iqr].min()
    
    outliers = np.append(v[v>=seventyfifth+1.5*iqr],v[v<=twentyfifth-1.5*iqr])
    
    #index median box_top box_bottom whisker_top whisker_bottom
    print("\\addplot+[boxplot prepared={"
        "draw position=",noise, 
        ",\n\tlower whisker=", lower_whisker, 
        ",\n\tlower quartile=", twentyfifth,
        ",\n\tmedian=", median,
        ",\n\tupper quartile=", seventyfifth, 
        ",\n\tupper whisker=", upper_whisker,
        "\n}]\ntable[row sep=\\\\,y index=0]{\n\t",
        '\\\\ '.join(str(i) for i in outliers)  ,
        "\\\n};",
        "\n\tcoordinates {};", sep="") 

v2 = np.array(np.mat(
'3.353483 3.429029 3.248354 3.701182 3.766603 8.637711 3.790586 3.586661  3.25941 3.260397 3.266955 3.157279 3.637177 4.720529 3.607241 5.141747 3.304324 4.176828  9.43966  3.74472 3.144122 4.314434 4.760472 3.245322 4.620148 3.420247  4.66521 5.064652 3.554964  3.40353 5.676936 3.329021 5.411193 3.467871  3.25708 3.539505 4.115666 3.381865 3.661999 4.192262 5.115789 3.502499 3.361452 3.920747 4.864098 6.415684 3.378859 3.344676 4.916818 3.371321 3.168973  3.62257 3.470234 3.389928 5.959497 3.723788 3.657174 4.491608 3.316763 3.459943  3.35802 3.219558 8.601731 3.267872 4.548442 9.264716 3.183346 3.480012 3.365831 3.272288 3.495957 3.389786 3.937944 6.017288  3.21402 3.933015 3.209189 3.396246 3.200412  4.42108 4.703574  3.27894 4.877642  3.90452 5.015378 3.527472 5.283496 3.133269 5.028169 3.638672 9.148105 3.488267 7.213787 3.572337 3.343015 8.962178 3.485171 9.828955 5.515183 3.317946'
));

for i in range(v2.shape[0]):
    format_data(v2[i,:], i+1)
    