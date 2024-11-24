export function fromHexString(hex){
    return Uint8Array.from(hex.match(/.{1,2}/g).map((byte) => parseInt(byte, 16)));
}

export function toHexString(bytes){
    return bytes.reduce((str, byte) => str + byte.toString(16).padStart(2, '0'), '');
}

export function toBase64(bytes) {
    return btoa(String.fromCharCode.apply(null, bytes));
}
  
export function fromBase64(base64) {
    return Uint8Array.from(atob(base64), c => c.charCodeAt(0))
}

export function toBase64URL(bytes) {
    const base64Encoded = toBase64(bytes);
    return base64Encoded.replace(/\+/g, '-').replace(/\//g, '_').replace(/=+$/, '');
}
  
export function fromBase64URL(base64) {
    const base64Encoded = base64.replace(/-/g, '+').replace(/_/g, '/');
    const padding = base64.length % 4 === 0 ? '' : '='.repeat(4 - (base64.length % 4));
    return fromBase64(base64Encoded + padding)
}
  