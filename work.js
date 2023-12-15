importScripts("gps-sdr-sim.js");

function sendProgress(percent){
    postMessage({type: 'progress', data: percent});
}
self.addEventListener('message', async (e)=> {
    console.log("worker run");
    const module = await SDRSim({
        exportProgress: sendProgress
    });
    module._staticMode(22.561092, 113.889375, 10, 0, 5.0, 8, 2600000, false);
    const file = module.FS.readFile('gpssim.bin', { encoding:'binary'});
    const blob = new Blob([file], { type: 'application/octet-stream' });
    const url = URL.createObjectURL(blob);
    postMessage({type: 'file', data: url});
})