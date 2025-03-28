% Ultralytics 🚀 AGPL-3.0 License - https://ultralytics.com/license

function [] = fcnAURA()
clc; close all; clear;
load('GE Synthetic - Vesuvius 1280x720 45FOV metadata.mat');

fig(1,1,4,3);  axis equal vis3d; set(gca,'zdir','reverse'); fcnview('skew')
fcnpaperplots(cam)
