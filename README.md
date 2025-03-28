<a href="https://www.ultralytics.com/"><img src="https://raw.githubusercontent.com/ultralytics/assets/main/logo/Ultralytics_Logotype_Original.svg" width="320" alt="Ultralytics logo"></a>

# MSVM (Minimum Separation Vector Mapping)

[![Ultralytics Actions](https://github.com/ultralytics/msvm/actions/workflows/format.yml/badge.svg)](https://github.com/ultralytics/msvm/actions/workflows/format.yml)
[![Ultralytics Discord](https://img.shields.io/discord/1089800235347353640?logo=discord&logoColor=white&label=Discord&color=blue)](https://discord.com/invite/ultralytics)
[![Ultralytics Forums](https://img.shields.io/discourse/users?server=https%3A%2F%2Fcommunity.ultralytics.com&logo=discourse&label=Forums&color=blue)](https://community.ultralytics.com/)
[![Ultralytics Reddit](https://img.shields.io/reddit/subreddit-subscribers/ultralytics?style=flat&logo=reddit&logoColor=white&label=Reddit&color=blue)](https://reddit.com/r/ultralytics)

Welcome to the Minimum Separation Vector Mapping (MSVM) project! This repository contains the implementation of an innovative [machine learning](https://www.ultralytics.com/glossary/machine-learning-ml) approach developed by Ultralytics for geospatial information fusion and video analytics. MSVM is designed to enhance situational awareness in intelligence, surveillance, and reconnaissance (ISR) applications, showcasing early work in advanced [computer vision](https://www.ultralytics.com/glossary/computer-vision-cv) techniques.

## 📑 Description

The MSVM technique, originally presented in our SPIE Defense + Security 2014 paper, utilizes sophisticated algorithms to map and analyze motion imagery specifically for ISR tasks. This method focuses on fusing geospatial data with video streams to provide deeper insights. For a comprehensive understanding, please consult the original publication:

Jocher, G., et al. "Minimum Separation Vector Mapping (MSVM)." Proc. SPIE 9089, Geospatial InfoFusion and Video Analytics IV; and Motion Imagery for ISR and Situational Awareness II, 90890A (2014). [DOI: 10.1117/12.2053833](https://doi.org/10.1117/12.2053833)

## 🔧 Requirements

To run the MSVM codebase, you need [MATLAB](https://www.mathworks.com/products/matlab.html) version 2018a or newer, along with specific toolboxes. Follow these setup steps:

1.  **Clone Common Functions:** Get the Ultralytics common MATLAB functions repository:
    ```shell
    git clone https://github.com/ultralytics/functions-matlab
    ```

2.  **Add to MATLAB Path:** Add the cloned repository directory to your MATLAB environment path using this command in MATLAB:
    ```matlab
    addpath(genpath('/path/to/functions-matlab'))
    ```
    Replace `/path/to/functions-matlab` with the actual path where you cloned the repository.

3.  **Install Required Toolboxes:** Ensure the following MATLAB toolboxes are installed:
    -   `Statistics and Machine Learning Toolbox`
    -   `Signal Processing Toolbox`

These toolboxes provide essential functions used by the MSVM algorithms.

## ▶️ Running the Code

To execute the MSVM estimators, open MATLAB, navigate to the project directory, and run the main script:

```matlab
runEstimators
```

This command will start the MSVM analysis process, generating output related to geospatial information fusion and video analytics based on the provided data.

Here is an example visualization of the expected results:
<img src="https://github.com/ultralytics/msvm/raw/main/results.jpg" alt="MSVM Results">

## 🤝 Contribute

We thrive on community contributions! If you're interested in fixing bugs, adding features, or improving documentation for MSVM or other [Ultralytics YOLO](https://docs.ultralytics.com/models/yolov8/) projects, your help is greatly appreciated. Please see our [Contributing Guide](https://docs.ultralytics.com/help/contributing/) for more details on how to get started.

We also value your feedback on Ultralytics products and encourage you to share your experiences by completing our brief [Survey](https://www.ultralytics.com/survey?utm_source=github&utm_medium=social&utm_campaign=Survey). A huge 🙏 thank you to all our contributors for supporting open-source AI! Explore more cutting-edge [AI solutions](https://docs.ultralytics.com/solutions/) and manage your projects with [Ultralytics HUB](https://www.ultralytics.com/hub).

[![Ultralytics open-source contributors](https://raw.githubusercontent.com/ultralytics/assets/main/im/image-contributors.png)](https://github.com/ultralytics/ultralytics/graphs/contributors)

## ©️ License

Ultralytics provides two licensing options to accommodate different use cases:

-   **AGPL-3.0 License**: This [OSI-approved](https://opensource.org/license/agpl-v3) open-source license is ideal for students, researchers, and enthusiasts keen on open collaboration and knowledge sharing. See the [LICENSE](https://github.com/ultralytics/msvm/blob/main/LICENSE) file for details.
-   **Enterprise License**: Designed for commercial applications, this license permits the seamless integration of Ultralytics software and AI models into commercial products and services, bypassing the open-source requirements of AGPL-3.0. If your project requires commercial licensing, please reach out through [Ultralytics Licensing](https://www.ultralytics.com/license).

## 📬 Contact Us

For bug reports, feature suggestions, or contributions related to the MSVM project, please visit [GitHub Issues](https://github.com/ultralytics/msvm/issues). For broader questions and discussions about Ultralytics projects and the AI community, join our active [Discord](https://discord.com/invite/ultralytics) server!

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
