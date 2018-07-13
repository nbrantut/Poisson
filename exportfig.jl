function findtextobj(x)
    return haskey(x,:set_fontsize)
end


function exportfig(fig, figname;
                   xsize=8.4,
                   ysize=6.5,
                   fontname="Times New Roman",
                   fontsize=9,
                   pad=.2,
                   w_pad=0.0,
                   h_pad=0.0)
    
    if fontname=="Times New Roman"
        matplotlib[:rc]("font", family="serif", serif="Times")
        matplotlib[:rc]("text.latex", preamble="\\usepackage{newtxmath}")
    elseif fontname=="Arial"
        matplotlib[:rc]("font", family="sans-serif")
        matplotlib[:rc]("text.latex", preamble="\\usepackage{newtxsf}")
    end

    figure(fig[:number])
    
    for item in fig[:findobj](findtextobj)
        item[:set_fontsize](fontsize)
        item[:set_fontname](fontname)
    end

    fig[:set_size_inches](xsize/2.54, ysize/2.54)

    

    fig[:canvas][:draw]()
        
    tight_layout(pad=pad, w_pad=w_pad, h_pad=h_pad)
    savefig(figname*".pdf")


end
