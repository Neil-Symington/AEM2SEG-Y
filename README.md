# AEM2SEG-Y
The AEM2SEG-Y package is a julia library that converts AEM conductivity line data into the SEG-Y format, which is commonly used for storing geophysical datasets such as reflection seismic or GPR. This allows AEM conductivity models to be imported visualised and ultimately interpreted within seismic software.

Feedback, suggestions or contributions will be gratefully accepted.

License
The content of this repository is licensed for use under the Apache 2.0 License.


Instructions for using this are here:

https://docs.google.com/document/d/1Elv19V3QclRzz9MhVqGacGvdHT4xERTr6etYUTFTpcA/edit?usp=sharing

This package is written in both python and Julia. Interestingly python is faster in the cases I tested as it seems to read text files much more quickly compared to Julia. Both languages rely on the segyio package.
Contacts

Neil Symington
neil.symington@ga.gov.au
