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

    tiledlayout(num,mum,'TileSpacing','compact','Padding','compact')
     

    % --- Sample ----------------------------------------------------------
    ax_q = nexttile(0*mum+1,[1 mum-1]);
    ax_p = nexttile(1*mum+1,[num-1 mum-1]);

    Gim.q = line(ax_q,Gim.sample_t,nan(size(Gim.sample_t)));
    xlim(ax_q,[chain.params.t_min chain.params.t_max]+0.1*(chain.params.t_max-chain.params.t_min)*[-1 +1])
    xline(ax_q,chain.params.t_min,'--','label','t_{min}','LabelOrientation','horizontal','LabelHorizontalAlignment','right')
    xline(ax_q,chain.params.t_max,'--','label','t_{max}','LabelOrientation','horizontal','LabelHorizontalAlignment','left')
    ylabel(ax_q,['nutrient density (',chain.params.units.conc,')'])

    Gim.p  = line(ax_p,Gim.sample_t,nan(size(Gim.sample_t)));
    Gim.c1 = line(ax_p,Gim.sample_t,nan(size(Gim.sample_t)),'linestyle','--');
    Gim.c2 = line(ax_p,Gim.sample_t,nan(size(Gim.sample_t)),'linestyle','--');
    xlim(ax_p,[chain.params.t_min chain.params.t_max]+0.1*(chain.params.t_max-chain.params.t_min)*[-1 +1])
    xline(ax_p,chain.params.t_min,'--','label','t_{min}','LabelOrientation','horizontal','LabelHorizontalAlignment','right')
    xline(ax_p,chain.params.t_max,'--','label','t_{max}','LabelOrientation','horizontal','LabelHorizontalAlignment','left')
    ylabel(ax_p,['cell density (',chain.params.units.conc,')'])

    xlabel(ax_p,['time (',chain.params.units.time,')'])

    title(ax_q,'MCMC sample')
    subtitle(ax_q,['(N=',num2str(chain.params.N),',K=',num2str(chain.params.K),')'])


    Gim.y = line(ax_p,repmat(chain.params.t,chain.params.K,1),nan(chain.params.N*chain.params.K,1),'marker','o','linestyle','none','color','c');
    Gim_z = line(ax_p,chain.params.t,chain.params.z,'marker','*','linestyle','-','color','k');

    legend(ax_p,[Gim.p,Gim.y,Gim_z],'trajectory','measurements','observations','location','northoutside','orientation','horizontal','box','off','autoupdate','off')

    y_lim_fl = true;


    % --- MCMC ------------------------------------------------------------
    ax_P = nexttile(1*mum,[2 1]);
    ax_h = nexttile(3*mum,[1 1]);
    
    ax_P.YAxisLocation = 'Right';
    ax_h.YAxisLocation = 'Right';

    ax_P.XLim = i_lim;
    ax_h.XLim = i_lim;
    
    title(ax_P,'MCMC chain')
    subtitle(ax_P,['(stride=',num2str(chain.stride),')'])

    xlabel(ax_h, 'MCMC iteration (i)')

    ax_P.YGrid = 'on';
  
    ylabel(ax_P,'logP_{post} (nat)' )
    ylabel(ax_h,'h (1)' )
 
    if ~isempty( chain.ledger )
        xline(ax_P,chain.ledger(end,1));
        xline(ax_h,chain.ledger(end,1));
    end

    Gim.P = line(ax_P,chain_i,chain.P(1,:),'marker','.');
    Gim.h = line(ax_h,chain_i,chain.h     ,'marker','.');


    % --- ground ----------------------------------------------------------
    
    if isfield( chain.params,'ground' )
        ground_x = get_x(chain.params.ground.g,Gim.sample_t,chain.params.t_min,chain.params.t_max);
        yline(ax_h,chain.params.ground.h,'-','color','g','linewidth',2)
        line(ax_q,Gim.sample_t,ground_x(:,1),'color','g','linewidth',2);
        line(ax_p,Gim.sample_t,ground_x(:,2),'color','g','linewidth',2);
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

drawnow

if y_lim_fl
    ylim(ax_q,"manual")
    ylim(ax_p,"manual")
end

end % visualize
