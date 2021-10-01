d=ded_dedalus_data_dir;


mm={'20','21','22','23','24','25','30','31','32','33','34','35'};
for j=1:length(mm);
  nm=['gc/f6/f/' mm{j}];
  a=load([d '/' nm '/profile.mat']);a=a.a;

  nn=nm;nn(nn=='/')='-';
  nn=['~/pics/mdns/' nn];
  
  figure(2*j-1);preprint([6 3],14);clf;colour_lines;
  h=plot(a.x-a.X,a.hu,a.x-a.X,a.hb);
  set(h,'linewidth',2);
  legend(h,{'u','b'},'location','SW');
  axis([-38 -2 0 1]);
  ylabel('$h$','interpreter','latex');
  xlabel('$x$','interpreter','latex');set(gca,'position',[0.11 0.2 0.88 0.76]);
  print('-depsc2',[nn '-h.eps']);
  
  figure(2*j);preprint([6 3],14);clf;colour_lines;
  h=plot(a.x-a.X,a.wu,a.x-a.X,a.wb);
  set(h,'linewidth',2);
  legend(h,{'u','b'},'location','SW');
  axis([-38 -2 0 0.6]);
  ylabel('$w$','interpreter','latex');
  xlabel('$x$','interpreter','latex');set(gca,'position',[0.11 0.2 0.88 0.76]);
  print('-depsc2',[nn '-w.eps']);
end


