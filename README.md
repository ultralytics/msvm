<img src="https://storage.googleapis.com/ultralytics/UltralyticsLogoName1000×676.png" width="200">

# Introduction 📘

Welcome to the official repository of the Minimum Separation Vector Mapping (MSVM) project, which hosts the innovative MSVM algorithm developed by the ultralytics team. This project marks a key contribution to the field of machine learning, specifically within the context of geospatial information fusion and video analytics. If you're interested in understanding or utilizing this advanced approach to machine learning, you've come to the right place.

# Description 📝

This repository is dedicated to the MSVM code, as introduced in the technical paper presented at the prestigious SPIE Defense + Security symposium in 2014:

- G. Jocher, et al., "Minimum Separation Vector Mapping (MSVM)," Proceedings of SPIE 9089, Geospatial InfoFusion and Video Analytics IV; and Motion Imagery for ISR and Situational Awareness II, 90890A (2014). [Link to Publication](http://dx.doi.org/10.1117/12.2053833)

The MSVM algorithm is an advanced machine learning approach designed to enhance analytics in geospatial information processing. The repository includes MATLAB code, which embodies this innovative algorithm, making it accessible for use and evaluation by the research community and industry professionals.

# Requirements ✅

Before you can run the MSVM code on your local machine, make sure you have the following requirements installed:

- [MATLAB](https://www.mathworks.com/products/matlab.html) version 2018a or newer.
- Download the common functions repository:

  ```bash
  $ git clone https://github.com/ultralytics/functions-matlab
  ```

  After cloning, add it to your MATLAB path:

  ```matlab
  >> addpath(genpath('/functions-matlab'))
  ```

- Additionally, ensure that you have the following MATLAB toolboxes installed:
  - `Statistics and Machine Learning Toolbox`
  - `Signal Processing Toolbox`

Please note that the code is meant to be compatible with MATLAB versions from 2018a onwards, and performance with earlier versions is not guaranteed.

# Running the Code 🚀

To execute the MSVM algorithm within MATLAB, navigate to the directory containing the code and run the following command:

```matlab
>> runEstimators  % This command initializes the MSVM process
```

Please be aware that running machine learning algorithms can be resource-intensive. Make sure your system meets the requirements mentioned above and has sufficient computational power to handle the processing.

The output of the MSVM algorithm can be visualized in MATLAB, providing insights into the effectiveness of the method and its application to geospatial analytics.

<img src="https://github.com/ultralytics/msvm/blob/master/results.jpg" alt="MSVM Results"> 

---

**Note:** This project is licensed under the AGPL-3.0 license. If you have any questions or if you would like to request permission for uses not covered by this license, please visit our website. 

For more information and updates, or if you are interested in further collaboration:

- Visit our website: [Ultralytics](http://www.ultralytics.com/contact)

Thank you for taking an interest in the MSVM project, and we look forward to seeing the innovative ways you can apply this algorithm!
