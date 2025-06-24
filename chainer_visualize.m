function Gim = chainer_visualize(Gim,chain)


y_lim_fl = false;


%% init
if isempty(Gim)
    
    num = 3;
    mum = 3;
    
    chain_i = double(chain.i(1)) + chain.stride*(0:chain.length-1)';
    i_lim = [max(chain_i(1),0.02*chain_i(end)) chain_i(end)+1];


    Gim.sample_t = linspace(chain.params.t_min,chain.params.t_max,10*chain.params.N)';
    
    figure(10)
    set(gcf,'windowstyle','docked')
    clf


    tiledlayout(num,mum)
     


    % --- Sample ----------------------------------------------------------
    ax_q = nexttile(0*mum+1,[1 mum-1]);
    ax_p = nexttile(1*mum+1,[num-1 mum-1]);


    Gim.q  = line(ax_q,Gim.sample_t,nan(size(Gim.sample_t)));
    xlim(ax_q,[chain.params.t_min chain.params.t_max]+0.1*(chain.params.t_max-chain.params.t_min)*[-1 +1])
    xline(ax_q,chain.params.t_min,'--','label','t_{min}','LabelOrientation','horizontal','LabelHorizontalAlignment','right')
    xline(ax_q,chain.params.t_max,'--','label','t_{max}','LabelOrientation','horizontal','LabelHorizontalAlignment','left')
    ylabel(ax_q,['nutrient density (',chain.params.units.conc,')'])

    Gim.p  = line(ax_p,Gim.sample_t,nan(size(Gim.sample_t)));
    Gim.c1 = line(ax_p,Gim.sample_t,nan(size(Gim.sample_t)),'linestyle','--','color',Gim.p.Color);
    Gim.c2 = line(ax_p,Gim.sample_t,nan(size(Gim.sample_t)),'linestyle','--','color',Gim.p.Color);
    xlim(ax_p,[chain.params.t_min chain.params.t_max]+0.1*(chain.params.t_max-chain.params.t_min)*[-1 +1])
    xline(ax_p,chain.params.t_min,'--','label','t_{min}','LabelOrientation','horizontal','LabelHorizontalAlignment','right')
    xline(ax_p,chain.params.t_max,'--','label','t_{max}','LabelOrientation','horizontal','LabelHorizontalAlignment','left')
    ylabel(ax_p,['cell density (',chain.params.units.conc,')'])


    xlabel(ax_p,['time (',chain.params.units.time,')'])


    title(ax_q,'MCMC sample')
    subtitle(ax_q,['(N=',num2str(chain.params.N),',K=',num2str(chain.params.K),')'])


    Gim.y = line(ax_p,repmat(chain.params.t,chain.params.K,1),nan(chain.params.N*chain.params.K,1),'marker','o','linestyle','none','color','c');
    Gim_z = line(ax_p,chain.params.t,chain.params.z(:,1),'marker','*','linestyle','-','color','k');


    if size(chain.params.z,2)==2
        line(ax_p,[chain.params.t';chain.params.t'],chain.params.z(:,1)' + chain.params.z(:,2)'.*[-1;+1],'marker','_','linestyle','-','color','k');
    end
    
    yline(0,':k')


    len = legend(ax_p,[Gim.p,Gim.y,Gim_z],'trajectory','measurements','observations','location','northoutside','orientation','horizontal','box','off');


    y_lim_fl = true;




    % --- MCMC ------------------------------------------------------------
    ax_P = nexttile(1*mum,[1 1]);
    ax_g = nexttile(2*mum,[1 1]);
    ax_h = nexttile(3*mum,[1 1]);


    ax_P.YAxisLocation = 'Right';
    ax_g.YAxisLocation = 'Right';
    ax_h.YAxisLocation = 'Right';

    ax_P.XLim = i_lim;
    ax_g.XLim = i_lim;
    ax_h.XLim = i_lim;
    
    title(ax_P,'MCMC chain')
    subtitle(ax_P,['(stride=',num2str(chain.stride),')'])


    xlabel(ax_h, 'MCMC iteration (i)')


    ax_P.YGrid = 'on';
    ax_g.YGrid = 'on';
    ax_h.YGrid = 'on';
  
    ax_g.YScale = 'log';

    ylabel(ax_P,'logP_{post} (nat)' )
    ylabel(ax_g,'g/\psi (1)' )
    ylabel(ax_h,'h (1)' )
 
    if ~isempty( chain.ledger )
        xline(ax_P,chain.ledger(end,1));
        xline(ax_g,chain.ledger(end,1));
        xline(ax_h,chain.ledger(end,1));
    end


    Gim.P = line(ax_P,chain_i,chain.P(1,:)                     ,'marker','.');
    Gim.g = line(ax_g,chain_i,chain.g./chain.params.g_prior_psi,'marker','.');
    Gim.h = line(ax_h,chain_i,chain.h                          ,'marker','.');




    % --- ground ----------------------------------------------------------
    
    if isfield( chain.params,'ground' )
        ground_x = get_x(chain.params.ground.g,Gim.sample_t,chain.params.t_min,chain.params.t_max);
        yline(ax_g,chain.params.ground.g./chain.params.g_prior_psi,'-','color','g','linewidth',2,'linestyle',':','label','ground','LabelHorizontalAlignment','left','LabelVerticalAlignment','Bottom')
        yline(ax_h,chain.params.ground.h                          ,'-','color','g','linewidth',2,'linestyle',':','label','ground','LabelHorizontalAlignment','left','LabelVerticalAlignment','Bottom')
        line(ax_q,Gim.sample_t,ground_x(:,1)                          ,'color','g','linewidth',2,'linestyle',':');
        line(ax_p,Gim.sample_t,ground_x(:,2)                          ,'color','g','linewidth',2,'linestyle',':','DisplayName','ground');
        len.AutoUpdate = 'off';
    end
    
end % init


sample_x = get_x(chain.sample.g,Gim.sample_t,chain.params.t_min,chain.params.t_max);

Gim.q.YData = sample_x(:,1);
Gim.p.YData = sample_x(:,2);

Gim.c1.YData = sample_x(:,2).*exp( sqrt(2/chain.sample.h)*erfinv( 2*0.05-1 ) );
Gim.c2.YData = sample_x(:,2).*exp( sqrt(2/chain.sample.h)*erfinv( 2*0.95-1 ) );


Gim.y.YData = chain.sample.y(:);


Gim.P.YData = chain.P(1,:);
Gim.h.YData = chain.h;
for i=1:length(chain.params.g_prior_psi)
    Gim.g(i).YData = chain.g(i,:)/chain.params.g_prior_psi(i);
end


drawnow


if y_lim_fl
    ylim(ax_q,"manual")
    ylim(ax_p,"manual")
end
