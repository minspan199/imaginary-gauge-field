classdef HologramHelperClass

    methods (Static)

        function E_Sample2 = supperposition(E_Sample)

            [M, N] = size(E_Sample);

            E_Sample2 = zeros(M, N);
            
            disp('Computing fields by Green function......')
            
            E_Sample3 = HologramHelperClass.supperposition2(E_Sample);
            
            disp('Assembling results......')
            
            for ii = 1:M

                for jj = 1:N    

                    E_Sample2 = E_Sample2 + E_Sample3{ii, jj};

                end

            end

        end
        
        function E_Sample3 = supperposition2(E_Sample)

            [M, N] = size(E_Sample);

            E_Sample2 = zeros(M, N);
            
            E_Sample3{M, N} = zeros(M, N);
            
            parfor ii = 1:M

                for jj = 1:N
                    
%                     if abs(E_Sample(ii, jj)) < 1e-6
%                         
%                         E_Sample3{ii, jj} = zeros(M, N);
%                         
%                         continue
%                         
%                     end                        

                    E_Sample3{ii, jj} = HologramHelperClass.Green(ii, jj, E_Sample);

                end

            end

        end

        function E_Sample = Green(coordImageX, coordImageY, E_Sample)

            global  k  N M z meshSize;
            meshSize = 2000E-9;
            lambda = 1.550e-6;
            k = 2 * pi / lambda;
            z = 1000e-6; % Distance of Image plane from the resonator plane
            dW = meshSize;
            [N, M] = size(E_Sample);
            amplitude = E_Sample(coordImageX, coordImageY); % inversely solved electric field at the sample plane

            SamplePlane = zeros(N, M); % electric field at sample plane.

            %     [coordSampleX, coordSampleY] = meshgrid((1:N) * dW, (1:M) * dW);

            for ii = 1:N

                for jj = 1:M

                    distance = (((coordImageX - ii) * dW)^2 + ((coordImageY - jj) * dW)^2 + z^2)^(-0.5);
                    SamplePlane(ii, jj) = amplitude * distance * exp(1i * k * distance);

                end

            end

            E_Sample = SamplePlane;
        end

        function normalizedIntensity = normalize(Original)

            [M, N] = size(Original);
            Emax = max(max(abs(Original)));

            for ii = 1:M

                for jj = 1:N

                    Original(ii, jj) = Original(ii, jj) / Emax;

                end

            end

            normalizedIntensity = Original;
        end
        
        
        function quantizedAmp = quantizeAmp(E, threshold)

            [M, N] = size(E);
            E = abs(E);
            
            for ii = 1:M

                for jj = 1:N

                    if(E(ii, jj) > 0.5)
                        E(ii, jj) = 1;
                    else
                        E(ii, jj) = 0;
                    end

                end

            end

            quantizedAmp = E;
        end
        
        
        function quantizedPhase = quantizePhase4(Array)

            [M, N] = size(Array);
            quantized = [-pi, -pi / 2, 0, pi / 2, pi]; % four level quantization
            % quantized = [-pi,-3*pi / 4 -pi / 2, -pi / 4, 0, pi / 4, pi / 2,3*pi / 4 pi]; % eight level quantization

            for ii = 1:M

                for jj = 1:N
                    [~, idx] = min(abs(quantized - Array(ii, jj)));
                    Array(ii, jj) = quantized(idx);

                end

            end

            quantizedPhase = Array;
        end
        
        function quantizedPhase = quantizePhase3(Array)

            [M, N] = size(Array);
            quantized = [-2*pi/3 0 2*pi/3]; % four level quantization
            % quantized = [-pi,-3*pi / 4 -pi / 2, -pi / 4, 0, pi / 4, pi / 2,3*pi / 4 pi]; % eight level quantization

            for ii = 1:M

                for jj = 1:N
                    [~, idx] = min(abs(quantized - Array(ii, jj)));
                    Array(ii, jj) = quantized(idx);

                end

            end

            quantizedPhase = Array;
        end
        
        function quantizedPhase = quantizePhase2(Array)

            [M, N] = size(Array);
            quantized = [0, pi]; % two level quantization
            % quantized = [-pi,-3*pi / 4 -pi / 2, -pi / 4, 0, pi / 4, pi / 2,3*pi / 4 pi]; % eight level quantization

            for ii = 1:M

                for jj = 1:N
                    [~, idx] = min(abs(quantized - Array(ii, jj)));
                    Array(ii, jj) = quantized(idx);

                end

            end

            quantizedPhase = Array;
        end

    end

end
