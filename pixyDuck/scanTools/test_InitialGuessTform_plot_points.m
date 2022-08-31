%%
figure(1); clf;
tiledlayout(1,2)


nexttile(1); hold on; axis equal;
title('rstageXY_fixed (g), stageXY_move (r)','Interpreter','none')
xlabel('x')
ylabel('y')

plot(rstageXY_fixed(:,1),rstageXY_fixed(:,2),'.','Color','g')

labelsFixed=arrayfun(@num2str,[1:size(rstageXY_fixed,1)]','UniformOutput',0);
text(rstageXY_fixed(:,1),rstageXY_fixed(:,2),labelsFixed)

plot(stageXY_move(:,1),stageXY_move(:,2),'.','Color','r')
labelsMove=arrayfun(@num2str,[1:size(stageXY_move,1)]','UniformOutput',0);
subset=[1,2,313];
text(stageXY_move(subset,1),stageXY_move(subset,2),labelsMove(subset))


%%
nexttile(2); hold on; axis equal;
title('withInitialGuessTform')
plot(rstageXY_fixed(:,1),rstageXY_fixed(:,2),'.','Color','g')
hold on
labelsFixed=arrayfun(@num2str,[1:size(rstageXY_fixed,1)]','UniformOutput',0);
text(rstageXY_fixed(:,1),rstageXY_fixed(:,2),labelsFixed)

plot(stageXY_move_WithInitialGuessTform(:,1),stageXY_move_WithInitialGuessTform(:,2),'.','Color','r')
labelsMove=arrayfun(@num2str,[1:size(stageXY_move_WithInitialGuessTform,1)]','UniformOutput',0);
subset=[1,2,313];
text(stageXY_move_WithInitialGuessTform(subset,1),stageXY_move_WithInitialGuessTform(subset,2),labelsMove(subset))

