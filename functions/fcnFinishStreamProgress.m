% Ultralytics 🚀 AGPL-3.0 License - https://ultralytics.com/license

function fcnFinishStreamProgress(h)
percentdone = get(h,'StreamingProgressPercentage');
if percentdone==100
    return
else
    while percentdone<100 %wait for stuff to render
        pause(.03)
        fprintf('streaming %.0f\n',percentdone)
        percentdone = get(h,'StreamingProgressPercentage');
    end
end

end

