function [J grad] = nnCostFunction(nn_params, ...
                                   input_layer_size, ...
                                   hidden_layer_size, ...
                                   num_labels, ...
                                   X, y, lambda)
%NNCOSTFUNCTION Implements the neural network cost function for a two layer
%neural network which performs classification
%   [J grad] = NNCOSTFUNCTON(nn_params, hidden_layer_size, num_labels, ...
%   X, y, lambda) computes the cost and gradient of the neural network. The
%   parameters for the neural network are "unrolled" into the vector
%   nn_params and need to be converted back into the weight matrices. 
% 
%   The returned parameter grad should be a "unrolled" vector of the
%   partial derivatives of the neural network.
%

% Reshape nn_params back into the parameters Theta1 and Theta2, the weight matrices
% for our 2 layer neural network
Theta1 = reshape(nn_params(1:hidden_layer_size * (input_layer_size + 1)), ...
                 hidden_layer_size, (input_layer_size + 1));

Theta2 = reshape(nn_params((1 + (hidden_layer_size * (input_layer_size + 1))):end), ...
                 num_labels, (hidden_layer_size + 1));

% Setup some useful variables
m = size(X, 1);
         
% You need to return the following variables correctly 
J = 0;
Theta1_grad = zeros(size(Theta1));
Theta2_grad = zeros(size(Theta2));

% ====================== YOUR CODE HERE ======================
% Instructions: You should complete the code by working through the
%               following parts.
%
% Part 1: Feedforward the neural network and return the cost in the
%         variable J. After implementing Part 1, you can verify that your
%         cost function computation is correct by verifying the cost
%         computed in ex4.m
%
% Part 2: Implement the backpropagation algorithm to compute the gradients
%         Theta1_grad and Theta2_grad. You should return the partial derivatives of
%         the cost function with respect to Theta1 and Theta2 in Theta1_grad and
%         Theta2_grad, respectively. After implementing Part 2, you can check
%         that your implementation is correct by running checkNNGradients
%
%         Note: The vector y passed into the function is a vector of labels
%               containing values from 1..K. You need to map this vector into a 
%               binary vector of 1's and 0's to be used with the neural network
%               cost function.
%
%         Hint: We recommend implementing backpropagation using a for-loop
%               over the training examples if you are implementing it for the 
%               first time.
%
% Part 3: Implement regularization with the cost function and gradients.
%
%         Hint: You can implement this around the code for
%               backpropagation. That is, you can compute the gradients for
%               the regularization separately and then add them to Theta1_grad
%               and Theta2_grad from Part 2.
% T1:[25x401];T2:[10x26]
[r1 c1] = size(Theta1);
[r2 c2] = size(Theta2);
DVec = zeros(1,hidden_layer_size * (input_layer_size + 1)+num_labels * (hidden_layer_size + 1));
tempT2 = Theta2(:,2:c2);
tempT1 = Theta1(:,2:c1);
%{
for i = 1: m
  t = zeros(1,num_labels);
  t(y(i)) = 1; % T [5000x10]
  g1 = X(i,:); %[1x400] G1[5000x400]
  a1 = [1 g1]; %[1x401] A1[5000x401]
  z2 = a1*Theta1'; %[1x25] Z2[5000x25]
  g2 = sigmoid(z2); %[1x25]
  a2 = [1 g2]; %[1x26] A2[5000X26]
  z3 = a2*Theta2'; %[1x10] Z3[5000x10]
  a3 = sigmoid(z3); %A3[5000x10]
  J = J + log(a3).*(-t)-log(1-a3).*(1-t);
  e3 = (a3 - t)'; %[10x1] E3[10x5000]
  e2 = tempT2'*e3.*g2'.*(1-g2'); %[25x1] E2[25x5000]
  e1 = tempT1'*e2.*g1'.*(1-g1'); %[400x1] E1[400x5000]
  d3 = e3*a2; %[10x26];
  d2 = e2*a1; %[25x401];
  DVec = DVec + [d2(:);d3(:)];
  fprintf(['%d\n'],i)
end
%}
Y = eye(num_labels)(y',:);
G1 = X;
A1 = [ones(m,1) G1];
Z2 = A1*Theta1';
G2 = sigmoid(Z2);
A2 = [ones(m,1) G2];
Z3 = A2*Theta2';
A3 = sigmoid(Z3);
E3 = (A3 - Y)';
E2 = tempT2'*E3.*G2'.*(1-G2');
E1 = tempT1'*E2.*G1'.*(1-G1');
D2 = E3*A2;
D1 = E2*A1;
DVec = [D1(:);D2(:)];
J = log(A3).*(-Y)-log(1-A3).*(1-Y); %[5000x10]



t1 = [zeros(r1,1) tempT1];
t2 = [zeros(r2,1) tempT2];
J = sum(sum(J))/m+lambda/2/m*(sum(sum(t1.^2))+sum(sum(t2.^2)));
D1 = reshape(DVec(1:hidden_layer_size * (input_layer_size + 1)), ...
                 hidden_layer_size, (input_layer_size + 1));
D2 = reshape(DVec((1 + (hidden_layer_size * (input_layer_size + 1))):end), ...
                 num_labels, (hidden_layer_size + 1));
Theta2_grad = D2/m + lambda*t2/m;
Theta1_grad = D1/m + lambda*t1/m;

% -------------------------------------------------------------

% =========================================================================

% Unroll gradients
grad = [Theta1_grad(:) ; Theta2_grad(:)];


end
