%Here alpha is computed similarly to the ChromaKey case. The main
%difference is given by the fact that the procedure is repeated for each
%cluster mean and covariance found in the window
function [bestF,bestB,bestAlpha] = computeAlpha(F_mi,F_sigma,B_mi,B_sigma,SIGMA,C,alpha)
loopN = 1;
previousL = inf;

for i=1:size(F_mi,2)
    F_miTemp = F_mi(:,i);
    F_sigmaTemp = F_sigma(:,:,i);
    for j=1:size(B_mi,2)
        B_miTemp = B_mi(:,j);
        B_sigmaTemp = B_sigma(:,:,j);
        while (loopN < 5)
            %matrix A
            first = inv(F_sigmaTemp)+eye(3)*(alpha^2)/SIGMA^2;
            second = eye(3)*alpha*(1-alpha)/SIGMA^2;
            third = eye(3)*alpha*(1-alpha)/SIGMA^2;
            fourth = inv(B_sigmaTemp)+eye(3)*(1-alpha^2)/SIGMA^2;
            A = [first, second; third, fourth];
            %vector b
            fifth = (F_sigmaTemp\F_miTemp)+C*alpha/SIGMA^2;
            sixth = (B_sigmaTemp\B_miTemp)+C*(1-alpha)/SIGMA^2;
            b = [fifth; sixth];
            %compute solution for system Ax = b
            x = A\b;
            %RGB values smaller tan 0 or bigger than 1 make no sense then
            %they are adjusted to stay in the range [0,1]
            F=max(0,min(1,x(1:3)));
            B=max(0,min(1,x(4:6)));
            alpha = max(0,min(1,dot(C-B,F-B)/dot(F-B,F-B)));            
            
            %compute the log likelihood for F anf B as described in the 
            %paper. L(alpha) is constant then is not considered for the 
            %computation of the log likelihood           
            LF=-((F-F_miTemp)'*inv(F_sigmaTemp)*(F-F_miTemp))/2;
            LB=-((B-B_miTemp)'*inv(B_sigmaTemp)*(B-B_miTemp))/2;
            L=LF+LB;
            
            %we save the likelihood with the biggest value
            if (previousL == inf || L>previousL)
                bestF = F;
                bestB = B;
                bestAlpha = alpha;
            end
                
            previousL = L;
            loopN = loopN +1;     
        end
        loopN = 1;
    end
end

