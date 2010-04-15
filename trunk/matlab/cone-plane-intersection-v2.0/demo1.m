function demo1
%
% Intersection between a cone and a plane.
% This source code is the implementation of the algorithms described in 
% Section 7.1.1, p.173 of the book "Robot Programming by Demonstration: A 
% Probabilistic Approach". 
%
% Author:	Sylvain Calinon, 2009
%			http://programming-by-demonstration.org
% 
% This program computes the intersection between a cone and a plane, and
% represents the intersection as a Gaussian Probability Density Function
% (PDF). The algorithm can be used to extract probabilistically information 
% concerning gazing or pointing direction. Indeed, by representing a visual 
% field as a cone and representing a table as a plane, the Gaussian PDF 
% can be used to compute the probability that one object on the table is 
% observed/pointed by the user. 
%
% This source code is given for free! However, I would be grateful if you refer 
% to the book (or corresponding article) in any academic publication that uses 
% this code or part of it. Here are the corresponding BibTex references: 
%
% @book{Calinon09book,
%   author="S. Calinon",
%   title="Robot Programming by Demonstration: A Probabilistic Approach",
%   publisher="EPFL/CRC Press",
%   year="2009",
%   note="EPFL Press ISBN 978-2-940222-31-5, CRC Press ISBN 978-1-4398-0867-2"
% }
% 
% @InProceedings{Calinon06roman,
%   author = "S. Calinon and A. Billard",
%   title = "Teaching a Humanoid Robot to Recognize and Reproduce Social Cues",
%   booktitle = "Proc. {IEEE} Intl Symposium on Robot and Human Interactive Communication ({Ro-Man})",
%   year = "2006",
%   month="September",
%   location="Hatfield, UK",
%   pages="346--351"
% }

%% Cone parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Cone vertex point
coneOrg = [0;0;2.0];
%Cone direction (axis of revolution), defined as a unity vector
%coneDir = [sin(-1.3)*cos(-1.1);-sin(-1.3)*cos(-1.1);sin(-1.1)];
coneDir = [0;0;-1];
%Cone half-angle
coneAngle = 0.2;
%Compute vectors (coneE1,coneE2) describing two orthogonal directions of 
%the cone basis.
coneE1 = cross(coneDir,[1;0;0]);
coneE2 = cross(coneDir,coneE1);

%% Plane parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Plane origin
planeOrg = [0;0;0];
%Plane first direction, defined as a unity vector
planeDir1 = [1;0;0];
%Plane second direction, defined as a unity vector
planeDir2 = [0;1;0];

%% Compute the cone-plane intersection represented as a Gaussian PDF
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[Mu, Sigma] = conePlaneIntersection(coneOrg, coneDir, coneAngle, planeOrg, planeDir1, planeDir2)

%% Subplot for the cone-plane intersection
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
subplot(1,2,1); hold on;
%Plot plane
planeSeg = [planeOrg-planeDir1-planeDir2 ...
  planeOrg+planeDir1-planeDir2 ...
  planeOrg+planeDir1+planeDir2 ...
  planeOrg-planeDir1+planeDir2 ...
  planeOrg-planeDir1-planeDir2];
patch(planeSeg(1,:),planeSeg(2,:),planeSeg(3,:),[.9 .9 .9]);
%Create segments for the intersection ellipse
nbSeg = 20;
t = linspace(-pi,pi,nbSeg)';
X = [cos(t) sin(t)] * real(sqrtm(Sigma)) + repmat(Mu',nbSeg,1);
for i=1:nbSeg
  interSeg(:,i) = planeOrg + X(i,1).*planeDir1 + X(i,2).*planeDir2; 
end
%Plot cone
for i=2:nbSeg
  patch([coneOrg(1) interSeg(1,i-1) interSeg(1,i)]', ...
    [coneOrg(2) interSeg(2,i-1) interSeg(2,i)]', ...
    [coneOrg(3) interSeg(3,i-1) interSeg(3,i)]', [1 0 0]);
end
%Plot cone-plane intersection
plot3(interSeg(1,:),interSeg(2,:),interSeg(3,:),'-','lineWidth',3,'color',[1 1 0]);
set(gcf,'Renderer','zbuffer'); lightangle(-45,30); lighting phong;
axis([-1 1 -1 1 0 2]); view(3); grid on;
xlabel('X','fontsize',16); ylabel('Y','fontsize',16); zlabel('Z','fontsize',16);

%% Subplot for the PDF representation of the intersection
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(1,2,2); hold on;
%Create a uniform grid
szgrid = 20;
[xx,yy] = meshgrid(linspace(-1,1,szgrid), linspace(-1,1,szgrid));
%Compute the probability density function for each point in the grid
zz = gaussPDF([xx(:) yy(:)]', Mu, Sigma);
zz = reshape(zz, size(xx));
%Plot the probability density function
surf(xx,yy,zz, 'LineStyle', '-'); 
axis([-1 1 -1 1 0 max(max(zz))]); colormap Gray; view(3); grid on;
xlabel('X','fontsize',16); ylabel('Y','fontsize',16); zlabel('P(X,Y)','fontsize',16);

pause;
close all;

 
