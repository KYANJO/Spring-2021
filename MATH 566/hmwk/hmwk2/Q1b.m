clear; clc;

A = [-2 11 ; -10 5]

plot_svd(A);


% defining function plot_svd

function plot_svd(A)

	t=0:0.001:(2*pi);
	% Create the unit circe on the pre-image plane.
	x1 = cos(t); y1 = sin(t);
	% Create the ellipse on the image plane.
	w = A*[x1; y1]; x2 = w(1,:); y2 = w(2,:);
	% Obtain the SVD of A
	[U, S, V] = svd(A);
	
	% Create the right and left singular vectors

	v1 = V(:,1); v2 = V(:,2);
	u1 = U(:,1); u2 = U(:,2);
	w1 = S(1,1)*u1; w2 = S(2,2)*u2;

	% Plot the unit circle and the right singular vectors on pre-image plane.
	figure(1);
	axis([-1.2 1.2 -1.2 1.2]);
	plot(x1, y1);
	hold on;
	line_plot(v1, 'b');
	hold on;
	line_plot(v2, 'r');
	grid on;
	title('Unit ball');
	
	% Plot the ellipse and the scaled left signualr vectors on image plane.
	figure(2);
	axis([-3 3 -3 3]);
	plot(x2, y2);
	hold on;
	line_plot(w1,'b');
	hold on;
	line_plot(w2, 'r');
	grid on;
	title('Unit ball image under A');
	end 

function line_plot(v, color)

	t = 0:0.01:1; x = v(1)*t; y = v(2)*t; plot(x, y, color);
end

