classdef fujisakiModel
% This class implements the Fujisaki model for f0 synthesis
% 
% Properties:
%   t (1,:) double vector
%       Time range for the model, must be non-empty and numeric.
%   F0b (1,1) double
%       Baseline F0, must be numeric and positive.
%   Ap (:,1) double vector
%       F0 phrase component, must be non-empty and numeric.
%   T0 (:,1) double vector
%       Time location of the phrase component, must be non-empty and numeric.
%   alpha (1,1) double
%       Decay rate for the phrase component, must be non-empty and numeric.
%   Aa (:,1) double vector
%       F0 accent component, must be non-empty and numeric.
%   beta (1,1) double
%       Decay rate for the accent component, must be non-empty and numeric.
%   gamma (1,1) double
%       Maximum value for the accent component, must be non-empty and numeric.
%   T1 (:,1) double vector
%       Onset time for the accent component, must be non-empty and numeric.
%   T2 (:,1) double vector
%       Offset time for the accent component, must be non-empty and numeric.
%   F0bContour (1,:) double vector
%       The baseline F0 contour, must be non-empty and numeric.
%   F0pContour (:,:) double matrix
%       The F0 phrase component contour, must be non-empty and numeric.
%   F0aContour (:,:) double matrix
%       The F0 accent component contour, must be non-empty and numeric.
%   F0Contour (:,:) double matrix
%       The F0 contour, must be non-empty and numeric.
% 
% Methods:
%   fujisakiModel
%       Constructor for the fujisakiModel class.
%       Input:
%           p: A struct containing the property values for the model.
%   assembleContour
%       Calculates the F0 baseline, phrase and accent components, and
%       stores the results in the F0bContour, F0pContour, F0aContour, and
%       F0Contour properties, respectively.
%   plotF0Contour
%       Plots the F0 contour.
% 
% Example:
%   p = struct('t', -1:.001:2.5, 'F0b', 70, 'Ap', [.2 .1], 'T0', [-.1 .5], ...
%       'alpha', 2, 'Aa', [.1 .05], 'beta', 20, 'gamma', .9, 'T1', [0 2.5], ...
%       'T2', [.2 4]);
%   obj = fujisakiModel(p);
%   obj.assembleContour();
%   obj.plotF0Contour(); 


    properties
        t (1,:) {mustBeNonempty,mustBeNumeric} = -1:.001:2.5
        F0b (1,1) {mustBeNumeric, mustBePositive} = 70
        Ap (:,1) {mustBeNonempty,mustBeNumeric} = [.2 .1]
        T0 (:,1) {mustBeNonempty,mustBeNumeric} = [-.1 .5]
        alpha (1,1) {mustBeNonempty,mustBeNumeric} = [2]
        Aa (:,1) {mustBeNonempty,mustBeNumeric} = [.1 .05];
        beta (1,1) {mustBeNonempty,mustBeNumeric}= 20;
        gamma (1,1) {mustBeNonempty,mustBeNumeric}= .9;
        T1 (:,1) {mustBeNonempty,mustBeNumeric}= [0 2.5];
        T2 (:,1) {mustBeNonempty,mustBeNumeric}= [.2 4];
        F0bContour (1,:) {mustBeNonempty,mustBeNumeric} = [0];
        F0pContour (:,:) {mustBeNonempty,mustBeNumeric} = [0];
        F0aContour (:,:) {mustBeNonempty,mustBeNumeric} = [0];
        F0Contour (:,:) {mustBeNonempty,mustBeNumeric} = [0];

    end

    methods
        function obj = fujisakiModel(p)
            % Construct an instance of class Fujisaki Model
            %   Detailed explanation goes here
            obj.t = p.t;
            obj.F0b = p.F0b;
            obj.Ap = p.Ap;
            obj.T0 = p.T0;
            obj.alpha = p.alpha;
            obj.Aa = p.Aa;
            obj.beta = p.beta;
            obj.gamma = p.gamma;
            obj.T1=p.T1;
            obj.T2 = p.T2;
        end

        function self = assembleContour(self)
            % METHOD1 Summary of this method goes here
            %   Detailed explanation goes here

            % Baseline F0
            self.F0bContour = log(ones(1,length(self.t)) .* self.F0b);

            % F0 phrase component
            self.F0pContour = (self.Ap) ...
                .* (self.alpha.^2) .* (self.t-self.T0)...
                .* exp(-self.alpha .* (self.t-self.T0))...
                .* double((self.t - self.T0) > 0);

            % F0 accent component
            GAOns = min(1-(1+self.beta.*(self.t - self.T1)).*exp(-self.beta.*(self.t - self.T1)),self.gamma)...
                .* double((self.t - self.T1) > 0);
            GAOff = min(1-(1+self.beta.*(self.t - self.T2)).*exp(-self.beta.*(self.t - self.T2)),self.gamma)...
                .* double((self.t - self.T2) > 0);
            self.F0aContour = (self.Aa).*(GAOns-GAOff);

            self.F0Contour = self.F0bContour+sum(self.F0pContour)+sum(self.F0aContour);
            % Plot
            % plot(self.t,exp(self.F0bContour),"LineWidth",2,"LineStyle",":","Color",[0 0 0]); hold on;
            % plot(self.t,exp(self.F0bContour+sum(self.F0pContour)),"LineWidth",2,"LineStyle","--","Color",[0 0 0])
            % plot(self.t,exp(self.F0bContour+sum(self.F0pContour)+sum(self.F0aContour)),"LineWidth",2,"LineStyle","-","Color",[0 0 0])
        end

        function f = plotF0Contour(self)
            f = figure;
            tiledlayout(5,1,"TileSpacing","compact")
            nexttile
            hold on
            plot(self.t,exp(self.F0bContour),"LineWidth",2,"LineStyle",":","Color",[0 0 0]);
            plot(self.t,exp(self.F0bContour+sum(self.F0pContour)),"LineWidth",2,"LineStyle","--","Color",[0 0 0])
            plot(self.t,exp(self.F0bContour+sum(self.F0pContour)+sum(self.F0aContour)),"LineWidth",2,"LineStyle","-","Color",[0 0 0])
            xlim([self.t(1) self.t(end)])
            ylabel("F0 [Hz]","Rotation",0,"HorizontalAlignment","right")
            grid on
            set(gca,"FontName","Assistant","FontSize",15)

            nexttile
            hold on
            for k = 1 : size(self.F0pContour,1)
                plot(self.t,exp(self.F0pContour(k,:)),"LineWidth",2,"LineStyle","--","Color",[0 0 0])
            end
            yline(1,"LineWidth",2,"LineStyle","-","Color",[0 0 0])
            xlim([self.t(1) self.t(end)])
            ylabel("Phrase Component","Rotation",0,"HorizontalAlignment","right")
            grid on
            set(gca,"FontName","Assistant","FontSize",15)

            nexttile
            hold on
            for k = 1 : size(self.F0pContour,1)
                plot([self.T0(k), self.T0(k)],[0,self.Ap(k)],"LineWidth",2,"LineStyle","-","Color",[0 0 0])
                markerSymb = "^";
                if self.Ap(k)<0
                    markerSymb = "v";
                end
                plot([self.T0(k)],[self.Ap(k)],markerSymb,"MarkerFaceColor",[0 0 0],"MarkerEdgeColor",[0 0 0])
            end
            yline(0,"LineWidth",2,"LineStyle","-","Color",[0 0 0])
            xlim([self.t(1) self.t(end)])
            ylabel("Phrase Command","Rotation",0,"HorizontalAlignment","right")
            grid on
            set(gca,"FontName","Assistant","FontSize",15)

            nexttile
            hold on
            for k = 1 : size(self.F0aContour,1)
                plot(self.t,exp(self.F0aContour(k,:)),"LineWidth",2,"LineStyle","-","Color",[0 0 0])
            end
            yline(1,"LineWidth",2,"LineStyle","-","Color",[0 0 0])
            xlim([self.t(1) self.t(end)])
            ylabel("Accent Component","Rotation",0,"HorizontalAlignment","right")
            grid on
            set(gca,"FontName","Assistant","FontSize",15)

            nexttile
            hold on
            for k = 1 : size(self.F0aContour,1)
                plot([self.T1(k), self.T1(k)],[0,self.Aa(k)],"LineWidth",2,"LineStyle","-","Color",[0 0 0])
                plot([self.T2(k), self.T2(k)],[0,self.Aa(k)],"LineWidth",2,"LineStyle","-","Color",[0 0 0])
                plot([self.T1(k), self.T2(k)],[self.Aa(k),self.Aa(k)],"LineWidth",2,"LineStyle","-","Color",[0 0 0])
            end
            yline(0,"LineWidth",2,"LineStyle","-","Color",[0 0 0])
            xlim([self.t(1) self.t(end)])
            ylabel("Accent Command","Rotation",0,"HorizontalAlignment","right")
            xlabel("Time [s]")
            grid on
            set(gca,"FontName","Assistant","FontSize",15)
        end
    end
end