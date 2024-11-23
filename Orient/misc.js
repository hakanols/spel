export function toBase64URL(bytes) {
    const base64Encoded = btoa(bytes);
    return base64Encoded.replace(/\+/g, '-').replace(/\//g, '_').replace(/=+$/, '');
}
  
export function fromBase64URL(str) {
    const base64Encoded = str.replace(/-/g, '+').replace(/_/g, '/');
    const padding = str.length % 4 === 0 ? '' : '='.repeat(4 - (str.length % 4));
    const base64WithPadding = base64Encoded + padding;
    return atob(base64WithPadding)
}

export function fromHexString(hexString){
    return Uint8Array.from(hexString.match(/.{1,2}/g).map((byte) => parseInt(byte, 16)));
}

export function toHexString(bytes){
    return bytes.reduce((str, byte) => str + byte.toString(16).padStart(2, '0'), '');
}
  