let isInit = false
init();

function init() {
    console.log("StartInit")
    isInit = true
    console.log("DoneInit")
};

export function getIsInit(){
    return isInit
}