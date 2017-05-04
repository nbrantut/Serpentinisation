function [I,check]=plane_line_intersect_vect(n_temp,V0_temp,P0,P1)
%plane_line_intersect computes the intersection of a plane and a segment(or
%a straight line)
% Inputs: 
%       n: normal vector of the Plane 
%       V0: any point that belongs to the Plane 
%       P0: end point 1 of the segment P0P1
%       P1:  end point 2 of the segment P0P1
%
%Outputs:
%      I    is the point of interection 
%     Check is an indicator:
%      0 => disjoint (no intersection)
%      1 => the plane intersects P0P1 in the unique point I
%      2 => the segment lies in the plane
%      3=>the intersection lies outside the segment P0P1
%
% Example:
% Determine the intersection of following the plane x+y+z+3=0 with the segment P0P1:
% The plane is represented by the normal vector n=[1 1 1]
% and an arbitrary point that lies on the plane, ex: V0=[1 1 -5]
% The segment is represented by the following two points
% P0=[-5 1 -1]
%P1=[1 2 3]   
% [I,check]=plane_line_intersect([1 1 1],[1 1 -5],[-5 1 -1],[1 2 3]);

%This function is written by :
%                             Nassim Khaled
%                             Wayne State University
%                             Research Assistant and Phd candidate
%If you have any comments or face any problems, please feel free to leave
%your comments and i will try to reply to you as fast as possible.
V0 = repmat(V0_temp,size(P1,1),1);
n  = repmat(n_temp,size(P1,1),1);
I=[0 0 0];
u = P1-P0;
w = P0 - V0;
D = dot(n,u,2);
N = -dot(n,w,2);
check=zeros(size(N));
check(D<10^-7 & N==0) = 2;
%compute the intersection parameter
sI = N./D;
I = P0+ repmat(sI,1,3).*u;
check(sI<0 | sI > 1) = 3; % The intersection is outside the segment so no intersection
check(sI>=0 & sI <= 1) = 1; % There is an intersection















