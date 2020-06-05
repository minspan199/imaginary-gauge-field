classdef HologramHelperClass

    methods (Static)

        function E_Sample2 = supperposition(E_Sample)

            [M, N] = size(E_Sample);

            E_Sample2 = zeros(M, N);

            for ii = 1:M

                for jj = 1:N

                    E_Sample2 = E_Sample2 + HologramHelperClass.Green(ii, jj, E_Sample);

                end

            end

        end

        function E_Sample = Green(coordImageX, coordImageY, E_Sample)

            global lambda k AmpImage N M dW z;
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

    end

end
