function [Mu, Sigma] = conePlaneIntersection(coneOrg, coneDir, coneAngle, planeOrg, planeDir1, planeDir2)
%
% This function computes the intersection between a cone and a plane, and
% represents the intersection as a Gaussian Probability Density Function
% (PDF). 
% This source code is the implementation of the algorithms described in 
% Section 7.1.1, p.173 of the book "Robot Programming by Demonstration: A 
% Probabilistic Approach". 
%
% Author:	Sylvain Calinon, 2009
%			http://programming-by-demonstration.org
% 
% The algorithm can be used to extract probabilistically information 
% concerning gazing or pointing direction. Indeed, by representing a visual 
% field as a cone and representing a table as a plane, the Gaussian PDF 
% can be used to compute the probability that one object on the table is 
% observed/pointed by the user.
%
% Inputs -----------------------------------------------------------------
%   o coneOrg:   3 x 1 vector representing the origin (vertex) of the cone  
%                in a Cartesian space.
%   o coneDir:   3 x 1 vector representing the direction (axis of 
%                revolution) of the cone in a Cartesian space.
%   o coneAngle: Value representing the half-cone angle in Radian.
%   o planeOrg:  3 x 1 vector representing the origin of the plane in a 
%                Cartesian space.
%   o planeDir1: 3 x 1 vector representing the first direction of the plane 
%                in a Cartesian space.
%   o planeDir2: 3 x 1 vector representing the second direction of the
%                plane in a Cartesian space.
% Outputs ----------------------------------------------------------------
%   o Mu:        2 x 1 vector representing the center of the Gaussian PDF 
%                on a 2D plane.
%   o Sigma:     2 x 2 array representing the covariance matrix of the 
%                Gaussian PDF on a 2D plane.
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

%% A point X on a cone is defined by the equation (in matricial form) 
%% (X-coneOrg)' * M * (X-coneOrg) = 0, 
%% with matrix M defined by:
M = coneDir * coneDir' - cos(coneAngle)^2.*diag(ones(3,1));

%% "Cone vertex - plane origin" directional vector 
conePlaneDir = planeOrg - coneOrg;

%% A point X on a plane is defined by equation
%% X = planeOrg + x1 * planeDir1 + x2 * planeDir2.
%% Combining the plane equation with the cone equation results in a
%% quadratic equation
%% c1*x1^2 + 2*c2*x1*x2 + c3*x2^2 + 2*c4*x1 + 2*c5*x2 + c6 = 0.
%% This equation can be represented as a matricial equation X'*C*X=0, 
%% with C the conic matrix.
%Elements of the conic matrix
c1 = planeDir1' * M * planeDir1;
c2 = planeDir1' * M * planeDir2;
c3 = planeDir2' * M * planeDir2;
c4 = conePlaneDir' * M * planeDir1;
c5 = conePlaneDir' * M * planeDir2;
c6 = conePlaneDir' * M * conePlaneDir;
%Conic matrix
C = [c1 c2 c4; c2 c3 c5; c4 c5 c6];
%Decomposition of the different parts of the conic matrix
C_R = C(1:2,1:2);
C_t = C(3,1:2)';
C_delta = C(3,3);

%% The intersection of a cone and a plane can form either an ellipse,
%% parabola or hyperbola. We are here only interested in the elliptical 
%% case.
%% The elliptical intersection conditions are:
if (det(C)~=0 && det(C_R)>0 && det(C)/(c1+c3)<0)==false
  disp('Error: the intersection is not an ellipse.');
  Mu = [0 0];
  Sigma = [1 0;0 1];
  return;
end

%% We determine a canonical form of the conic by transforming the
%% original conic through a rotation R and a translation t, i.e. by
%% applying an homogenous transformation H.
%% To find R, we diagonalize C_R by computing its eigenvalues, i.e.
%% C_R = R*Lambda*R'.
[R,Lambda] = eig(C_R);

%% The translation t centering the ellipse to the canonical form is 
%% determined by: 
t = -R * inv(Lambda) * R' * C_t;

%% With rotation R and translation t, we can then define the homogenous 
%% transformation matrix H
H = [R t; 0 0 1];

%% The canonical conic matrix C_c is then defined by
%% C_c = H' * C * H.
C_c = H' * C * H;

%% C_c is a diagonal matrix with elements (C_c1,C_c2,C_c3), which defines 
%% the canonical quadratic equation:
%% C_c1 * x_c1^2 + C_c2 * x_c2^2 + C_c3 = 0,
%% This quadratic equation can be re-written as an ellips equation:
%% x_c1^2 / a^2 + x_c2^2 / b^2 = 1,
%% with parameters:
a = sqrt(-C_c(3,3)/C_c(1,1));
b = sqrt(-C_c(3,3)/C_c(2,2));

%% Similarly, the ellipse can aso be represented as a diagonal 
%% covariance matrix Sigma_c: 
Sigma_c = [a^2 0; 0 b^2];

%% Finally, the 2D Gaussian distribution (Mu,Sigma) is computed using 
%% the translation vector t and the canonical covariance matrix Sigma_c 
%% transformed through rotation R.   
Mu = t;
Sigma = R * Sigma_c * R';

  
