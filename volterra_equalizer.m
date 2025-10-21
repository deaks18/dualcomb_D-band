function [x_hat, w, err] = volterra_equalizer(x, y, M1, M2, M3, mu1, mu2, mu3, order)
% VOLTERRA_EQUALIZER_CUSTOM - Volterra equalizer with different memory lengths and step sizes per order
%
% Inputs:
%   x     - Original transmitted signal
%   y     - Received distorted signal
%   M1    - Linear memory length
%   M2    - Quadratic memory length
%   M3    - Cubic memory length
%   mu1   - Step size for linear terms
%   mu2   - Step size for quadratic terms
%   mu3   - Step size for cubic terms
%   order - Maximum nonlinear order (1, 2, or 3)
%
% Outputs:
%   x_hat - Estimated output signal
%   w     - Final weight vector (all orders)
%   err   - Error signal

    N = length(x);
    x_hat = zeros(N, 1);
    err = zeros(N, 1);

    % Number of terms per order
    n1 = M1;
    n2 = (order >= 2) * (M2 * (M2 + 1) / 2);
    n3 = 0;
    if order == 3
        n3 = nchoosek(M3 + 2, 3);  % symmetric cubic terms
    end
    total_terms = n1 + n2 + n3;

    % Initialize weights
    w = zeros(total_terms, 1);

    % Index mapping
    idx1 = 1;
    idx2 = idx1 + n1;
    idx3 = idx2 + n2;

    for n = max([M1, M2, M3]) : N
        phi = zeros(total_terms, 1);

        %% Linear terms
        if M1 > 0
            buf1 = y(n:-1:n - M1 + 1);
            phi(idx1 : idx1 + n1 - 1) = buf1;
        end

        %% Quadratic terms (i ≤ j)
        if order >= 2 && M2 > 0
            buf2 = y(n:-1:n - M2 + 1);
            idx = 0;
            for i = 1:M2
                for j = i:M2
                    idx_phi = idx2 + idx;
                    phi(idx_phi) = buf2(i) * buf2(j);
                    idx = idx + 1;
                end
            end
        end

        %% Cubic terms (i ≤ j ≤ k)
        if order == 3 && M3 > 0
            buf3 = y(n:-1:n - M3 + 1);
            idx = 0;
            for i = 1:M3
                for j = i:M3
                    for k = j:M3
                        idx_phi = idx3 + idx;
                        phi(idx_phi) = buf3(i) * buf3(j) * buf3(k);
                        idx = idx + 1;
                    end
                end
            end
        end

        % Output estimate
        x_hat(n) = w' * phi;

        % Error
        err(n) = x(n) - x_hat(n);

        % Update weights by order
        if M1 > 0
            w(idx1 : idx1 + n1 - 1) = w(idx1 : idx1 + n1 - 1) + mu1 * err(n) * phi(idx1 : idx1 + n1 - 1);
        end
        if order >= 2 && M2 > 0
            w(idx2 : idx2 + n2 - 1) = w(idx2 : idx2 + n2 - 1) + mu2 * err(n) * phi(idx2 : idx2 + n2 - 1);
        end
        if order == 3 && M3 > 0
            w(idx3 : idx3 + n3 - 1) = w(idx3 : idx3 + n3 - 1) + mu3 * err(n) * phi(idx3 : idx3 + n3 - 1);
        end
    end
end
