<br>
<a href="https://www.ultralytics.com/" target="_blank"><img src="https://raw.githubusercontent.com/ultralytics/assets/main/logo/Ultralytics_Logotype_Original.svg" width="320" alt="Ultralytics logo"></a>

# 🌟 Introduction

Welcome to the MSVM (Minimum Separation Vector Mapping) project! This repository hosts the implementation of an innovative machine learning approach for geospatial information fusion and video analytics, designed to enhance situational awareness in intelligence, surveillance, and reconnaissance (ISR) applications.

[![Ultralytics Actions](https://github.com/ultralytics/msvm/actions/workflows/format.yml/badge.svg)](https://github.com/ultralytics/msvm/actions/workflows/format.yml) <a href="https://discord.com/invite/ultralytics"><img alt="Discord" src="https://img.shields.io/discord/1089800235347353640?logo=discord&logoColor=white&label=Discord&color=blue"></a> <a href="https://community.ultralytics.com/"><img alt="Ultralytics Forums" src="https://img.shields.io/discourse/users?server=https%3A%2F%2Fcommunity.ultralytics.com&logo=discourse&label=Forums&color=blue"></a> <a href="https://reddit.com/r/ultralytics"><img alt="Ultralytics Reddit" src="https://img.shields.io/reddit/subreddit-subscribers/ultralytics?style=flat&logo=reddit&logoColor=white&label=Reddit&color=blue"></a>

# 📑 Description

The MSVM technique, introduced in our SPIE Defense + Security 2014 paper, employs advanced algorithms to map and analyze motion imagery for ISR tasks. For more details, please refer to our publication:

Jocher, G., et al. "Minimum Separation Vector Mapping (MSVM)." Proc. SPIE 9089, Geospatial InfoFusion and Video Analytics IV; and Motion Imagery for ISR and Situational Awareness II, 90890A (2014). [DOI: 10.1117/12.2053833](http://dx.doi.org/10.1117/12.2053833)

# 🔧 Requirements

Before diving into the MSVM codebase, ensure that you have [MATLAB](https://www.mathworks.com/products/matlab.html) version 2018a or newer installed, alongside the necessary toolboxes for smooth operation. Follow these steps to get started:

1. Clone the common functions repository provided by Ultralytics:

```shell
$ git clone https://github.com/ultralytics/functions-matlab
```

2. Add the cloned repository to your MATLAB path by executing the following command in MATLAB:

```matlab
>> addpath(genpath('/path/to/functions-matlab'))
```

Make sure to have the following MATLAB toolboxes installed and ready:

- `Statistics and Machine Learning Toolbox`
- `Signal Processing Toolbox`

These toolboxes are essential for executing the ML algorithms within the MSVM framework.

# ▶️ Running the Code

To execute the MSVM estimators, simply navigate to the MATLAB environment and run the following command:

```matlab
>> runEstimators
```

This will initiate the MSVM analysis and provide you with the output related to geospatial information fusion and video analytics.

Below is a snapshot of the expected results: <img src="https://github.com/ultralytics/msvm/blob/main/results.jpg" alt="MSVM Results">

# 🤝 Contribute

We welcome contributions from the community! Whether you're fixing bugs, adding new features, or improving documentation, your input is invaluable. Take a look at our [Contributing Guide](https://docs.ultralytics.com/help/contributing/) to get started. Also, we'd love to hear about your experience with Ultralytics products. Please consider filling out our [Survey](https://www.ultralytics.com/survey?utm_source=github&utm_medium=social&utm_campaign=Survey). A huge 🙏 and thank you to all of our contributors!

<!-- Ultralytics contributors -->

<a href="https://github.com/ultralytics/yolov5/graphs/contributors">
<img width="100%" src="https://github.com/ultralytics/assets/raw/main/im/image-contributors.png" alt="Ultralytics open-source contributors"></a>

# ©️ License

Ultralytics is excited to offer two different licensing options to meet your needs:

- **AGPL-3.0 License**: Perfect for students and hobbyists, this [OSI-approved](https://opensource.org/license) open-source license encourages collaborative learning and knowledge sharing. Please refer to the [LICENSE](https://github.com/ultralytics/ultralytics/blob/main/LICENSE) file for detailed terms.
- **Enterprise License**: Ideal for commercial use, this license allows for the integration of Ultralytics software and AI models into commercial products without the open-source requirements of AGPL-3.0. For use cases that involve commercial applications, please contact us via [Ultralytics Licensing](https://www.ultralytics.com/license).

# 📬 Contact Us

For bug reports, feature requests, and contributions, head to [GitHub Issues](https://github.com/ultralytics/velocity/issues). For questions and discussions about this project and other Ultralytics endeavors, join us on [Discord](https://discord.com/invite/ultralytics)!

<br>
<div align="center">
  <a href="https://github.com/ultralytics"><img src="https://github.com/ultralytics/assets/raw/main/social/logo-social-github.png" width="3%" alt="Ultralytics GitHub"></a>
  <img src="https://github.com/ultralytics/assets/raw/main/social/logo-transparent.png" width="3%" alt="space">
  <a href="https://www.linkedin.com/company/ultralytics/"><img src="https://github.com/ultralytics/assets/raw/main/social/logo-social-linkedin.png" width="3%" alt="Ultralytics LinkedIn"></a>
  <img src="https://github.com/ultralytics/assets/raw/main/social/logo-transparent.png" width="3%" alt="space">
  <a href="https://twitter.com/ultralytics"><img src="https://github.com/ultralytics/assets/raw/main/social/logo-social-twitter.png" width="3%" alt="Ultralytics Twitter"></a>
  <img src="https://github.com/ultralytics/assets/raw/main/social/logo-transparent.png" width="3%" alt="space">
  <a href="https://youtube.com/ultralytics"><img src="https://github.com/ultralytics/assets/raw/main/social/logo-social-youtube.png" width="3%" alt="Ultralytics YouTube"></a>
  <img src="https://github.com/ultralytics/assets/raw/main/social/logo-transparent.png" width="3%" alt="space">
  <a href="https://www.tiktok.com/@ultralytics"><img src="https://github.com/ultralytics/assets/raw/main/social/logo-social-tiktok.png" width="3%" alt="Ultralytics TikTok"></a>
  <img src="https://github.com/ultralytics/assets/raw/main/social/logo-transparent.png" width="3%" alt="space">
  <a href="https://ultralytics.com/bilibili"><img src="https://github.com/ultralytics/assets/raw/main/social/logo-social-bilibili.png" width="3%" alt="Ultralytics BiliBili"></a>
  <img src="https://github.com/ultralytics/assets/raw/main/social/logo-transparent.png" width="3%" alt="space">
  <a href="https://discord.com/invite/ultralytics"><img src="https://github.com/ultralytics/assets/raw/main/social/logo-social-discord.png" width="3%" alt="Ultralytics Discord"></a>
</div>
