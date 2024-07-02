function [u_field,v_field] = DataWrapper()
    clear all;
    close all;
    format compact
    format long

    global U = [];
    global V = [];

    load u_2014.mat;   
    load v_2014.mat;  % These two datasets contain the velocity data in the whole year 2014.

    [u_field,v_field]=DataPreprocessing(u_2014,v_2014);
    U=u_field;
    V=v_field;