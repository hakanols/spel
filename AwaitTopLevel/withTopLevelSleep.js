let isInit = false
await init();

async function init() {
    console.log("StartInit")
    await sleep(500)
    isInit = true
    console.log("DoneInit")
};

export function getIsInit(){
    return isInit
}

function sleep(ms) {
    return new Promise(resolve => setTimeout(resolve, ms));
}