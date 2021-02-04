# AEM2SEG-Y
The AEM2SEG-Y package is a python library that converts AEM conductivity line data into the SEG-Y format, whcih is commonly used within for storing geophysical dataset such as reflection seismic or GOR. This allows AEM conductivity models to be imported visualised and ultimately interpreted within seismic software.

Feedback, suggestions or contributions will be gratefully accepted.

License
The content of this repository is licensed for use under the Apache 2.0 License.


Instructions for using this are here:

https://docs.google.com/document/d/1Elv19V3QclRzz9MhVqGacGvdHT4xERTr6etYUTFTpcA/edit?usp=sharing

Known bugs include the software throwing errors if the AEM line is too long (e.g. AusAEM survey). This is inherent segy and hence the lines will need to be chunked before conversion.

Contacts

Neil Symington
neil.symington@ga.gov.au
