function draw_box(x, y, side, color, w)

A = [x - side/2; y - side/2];
B = [x - side/2; y + side/2];
C = [x + side/2; y - side/2];
D = [x + side/2; y + side/2];

hold on;
plot([A(1) B(1)], [A(2) B(2)], color, 'LineWidth', w)
plot([D(1) B(1)], [D(2) B(2)], color, 'LineWidth', w)
plot([A(1) C(1)], [A(2) C(2)], color, 'LineWidth', w)
plot([C(1) D(1)], [C(2) D(2)], color, 'LineWidth', w)