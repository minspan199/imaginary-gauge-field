
vertical_bar12by12 = [0	0	0	0	0	0	0	0	0	0	0	0
0	0	0	0	0	0	0	0	0	0	0	0
0	0	0	0	0	0	0	0	0	0	0	0
0	0	0	0	0	0	0	0	0	0	0	0
0	0	0	0	0	0	1	0	0	0	0	0
0	0	0	0	0	0	1	0	0	0	0	0
0	0	0	0	0	0	1	0	0	0	0	0
0	0	0	0	0	0	1	0	0	0	0	0
0	0	0	0	0	0	0	0	0	0	0	0
0	0	0	0	0	0	0	0	0	0	0	0
0	0	0	0	0	0	0	0	0	0	0	0
0	0	0	0	0	0	0	0	0	0	0	0
];

vertical_bar10by10 = [0	0	0	0	0	0	0	0	0	0
0	0	0	0	0	0	0	0	0	0
0	0	0	0	0	0	0	0	0	0
0	0	0	0	1	0	0	0	0	0
0	0	0	0	1	0	0	0	0	0
0	0	0	0	1	0	0	0	0	0
0	0	0	0	1	0	0	0	0	0
0	0	0	0	0	0	0	0	0	0
0	0	0	0	0	0	0	0	0	0
0	0	0	0	0	0	0	0	0	0
];

vertical_bar8by8 = [0	0	0	0	0	0	0	0
0	0	0	1	0	0	0	0
0	0	0	1	0	0	0	0
0	0	0	1	0	0	0	0
0	0	0	1	0	0	0	0
0	0	0	1	0	0	0	0
0	0	0	1	0	0	0	0
0	0	0	0	0	0	0	0

];

circle1 = [1	1	0	0	0	0	0	0	0	0
1	1	1	1	1	1	1	1	1	0
1	1	0	0	0	0	0	0	1	0
0	1	0	0	0	0	0	0	1	0
0	1	0	0	0	0	0	0	1	0
0	1	0	0	0	0	0	0	1	0
0	1	0	0	0	0	0	0	1	0
0	1	0	0	0	0	0	0	1	0
0	1	1	1	1	1	1	1	1	0
0	0	0	0	0	0	0	0	0	0
];

horizontal_bar10by10 = [0	0	0	0	0	0	0	0	0	0
0	0	0	0	0	0	0	0	0	0
0	0	0	0	0	0	0	0	0	0
0	0	0	0	0	0	0	0	0	0
0	0	0	0	0	0	0	0	0	0
0	0	0	1	1	1	1	0	0	0
0	0	0	0	0	0	0	0	0	0
0	0	0	0	0	0	0	0	0	0
0	0	0	0	0	0	0	0	0	0
0	0	0	0	0	0	0	0	0	0
];

figure;
imagesc(horizontal_bar10by10, [0, 1]);
set(gcf, 'Position', [00, 00, 350, 300])
set(gca, 'FontSize', 10)% Font Size
colormap gray
colorbar
title('Desired Far Field Image');

imwrite(horizontal_bar10by10, 'horizontal_bar_10by10.png');

IntensityNear = ones(8, 8);
save('bar_8by8_raw.mat','IntensityNear');



%% FarField


% clear all

N = 20;


IntensityNear = [1	1	1	1	1	1	1	1	1	1
1	0	0	0	0	0	0	0	0	1
1	0	0	0	0	0	0	0	0	1
1	0	0	0	0	0	0	0	0	1
1	0	0	0	0	0	0	0	0	1
1	0	0	0	0	0	0	0	0	1
1	0	0	0	0	0	0	0	0	1
1	0	0	0	0	0	0	0	0	1
1	0	0	0	0	0	0	0	0	1
1	1	1	1	1	1	1	1	1	1
];

IntensityNear = [1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1
1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1
1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1
1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1
1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1
1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1
1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1
1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1
1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1
1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1
1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1
1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1
1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1
1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1
1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1
1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1
1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1
1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1
1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1
1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1
];

IntensityNear = [0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
0	0	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	0	0
0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0
0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0
0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0
0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0
0	0	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	0	0
0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
0	0	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	0	0
0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0
0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0
0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0
0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0
0	0	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	0	0
0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
];


save('imageNear_grids.mat', 'IntensityNear');

figure
imagesc(IntensityNear);

IntensityNear(N, N) = 0;
IntensityNear(1:2:N - 1, :) = 1;
IntensityNear(:, 1:2:N - 1) = 1;
save('imageNear_meshgrid.mat', 'IntensityNear');

IntensityNear = ones(N, N);
save('imageNear_uniform.mat', 'IntensityNear');

%% A letter map
A = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
    0 0 0 0 0 0 0 0 1 1 1 1 0 0 0 0 0 0 0 0
    0 0 0 0 0 0 0 1 0 0 0 0 1 0 0 0 0 0 0 0
    0 0 0 0 0 0 1 0 0 0 0 0 0 1 0 0 0 0 0 0
    0 0 0 0 0 1 0 0 0 0 0 0 0 0 1 0 0 0 0 0
    0 0 0 0 0 1 0 0 0 0 0 0 0 0 1 0 0 0 0 0
    0 0 0 0 0 1 0 0 0 0 0 0 0 0 1 0 0 0 0 0
    0 0 0 0 0 1 0 0 0 0 0 0 0 0 1 0 0 0 0 0
    0 0 0 0 0 1 0 0 0 0 0 0 0 0 1 0 0 0 0 0
    0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0
    0 0 0 0 0 1 0 0 0 0 0 0 0 0 1 0 0 0 0 0
    0 0 0 0 0 1 0 0 0 0 0 0 0 0 1 0 0 0 0 0
    0 0 0 0 0 1 0 0 0 0 0 0 0 0 1 0 0 0 0 0
    0 0 0 0 0 1 0 0 0 0 0 0 0 0 1 0 0 0 0 0
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
    ];
imwrite(A, 'LetterA_thin.png')

%% B letter map
B = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
    0 0 0 0 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0
    0 0 0 0 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0
    0 0 0 0 1 1 0 0 0 0 0 0 0 1 1 0 0 0 0 0
    0 0 0 0 1 1 0 0 0 0 0 0 0 1 1 0 0 0 0 0
    0 0 0 0 1 1 0 0 0 0 0 0 0 1 1 0 0 0 0 0
    0 0 0 0 1 1 0 0 0 0 0 0 0 1 1 0 0 0 0 0
    0 0 0 0 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0
    0 0 0 0 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0
    0 0 0 0 1 1 0 0 0 0 0 0 0 1 1 0 0 0 0 0
    0 0 0 0 1 1 0 0 0 0 0 0 0 1 1 0 0 0 0 0
    0 0 0 0 1 1 0 0 0 0 0 0 0 1 1 0 0 0 0 0
    0 0 0 0 1 1 0 0 0 0 0 0 0 1 1 0 0 0 0 0
    0 0 0 0 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0
    0 0 0 0 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
    ];

%% C letter map
C = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
    0 0 0 0 0 0 1 1 1 1 1 1 1 0 0 0 0 0 0 0
    0 0 0 0 0 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0
    0 0 0 0 1 1 0 0 0 0 0 0 0 1 1 0 0 0 0 0
    0 0 0 1 1 0 0 0 0 0 0 0 0 1 1 0 0 0 0 0
    0 0 0 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
    0 0 0 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
    0 0 0 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
    0 0 0 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
    0 0 0 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
    0 0 0 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
    0 0 0 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0
    0 0 0 0 1 1 0 0 0 0 0 0 0 1 1 0 0 0 0 0
    0 0 0 0 0 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0
    0 0 0 0 0 0 1 1 1 1 1 1 1 0 0 0 0 0 0 0
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
    ];

%% 45 degree bar

bar45 = [0	0	0	0	0	0	0	0	0	0
0	0	0	0	0	0	0	0	0	0
0	0	0	0	0	0	0	0	0	0
0	0	0	0	0	0	1	0	0	0
0	0	0	0	0	1	0	0	0	0
0	0	0	0	1	0	0	0	0	0
0	0	0	1	0	0	0	0	0	0
0	0	0	0	0	0	0	0	0	0
0	0	0	0	0	0	0	0	0	0
0	0	0	0	0	0	0	0	0	0
];

cross_45 = [0	0	0	0	0	0	0	0	0	0
0	0	0	0	0	0	0	0	0	0
0	0	1	0	0	0	0	1	0	0
0	0	0	1	0	0	1	0	0	0
0	0	0	0	1	1	0	0	0	0
0	0	0	0	1	1	0	0	0	0
0	0	0	1	0	0	1	0	0	0
0	0	1	0	0	0	0	1	0	0
0	0	0	0	0	0	0	0	0	0
0	0	0	0	0	0	0	0	0	0
];
imwrite(cross_45, 'cross_45.png')


cross_0 = [0	0	0	0	0	0	0	0	0	0
0	0	0	0	0	0	0	0	0	0
0	0	0	0	0	0	0	0	0	0
0	0	0	0	0	1	0	0	0	0
0	0	0	0	0	1	0	0	0	0
0	0	0	1	1	1	1	1	0	0
0	0	0	0	0	1	0	0	0	0
0	0	0	0	0	1	0	0	0	0
0	0	0	0	0	0	0	0	0	0
0	0	0	0	0	0	0	0	0	0

];
imwrite(cross_0, 'cross_0.png')

