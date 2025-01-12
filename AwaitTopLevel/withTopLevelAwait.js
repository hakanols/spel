await init();

let isInit = false

async function init() {
    console.log("StartInit")
    await sleep(2000)
    isInit = true
    console.log("DoneInit")
};

export function getIsInit(){
    return isInit
}

function sleep(ms) {
    return new Promise(resolve => setTimeout(resolve, ms));
}