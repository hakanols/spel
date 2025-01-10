import nacl from './nacl-fast-es.js';
import * as misc from './misc.js';

export function getStaticString() {
    const seed = (new TextEncoder()).encode("RandomText");
    return misc.toBase64URL(nacl.hash(seed).slice(0,4))
}

export function getRandomString() {
    return misc.toBase64URL(nacl.randomBytes(4))
}